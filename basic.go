package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"math/rand"
	"os"
	"strconv"
	"sync"
	"time"
)

// All SI units.  The genesis of every value is in SI units
// All densities are number densities m^-3
// Unit conversion only occur after computation
// Global constants
const (
	AvogadroNumber = 6.02214076e23     // Avogadro's number
	AMUtoKG        = 1.66053906660e-27 // Conversion factor from AMU to kg
	eVtoJ          = 1.602e-19         // Conversion from eV to J
	JtoeV          = 1 / eVtoJ         // Convertion from J to eV
	DeuteriumMass  = 2.014 * AMUtoKG   // masses are in kg
	TitaniumMass   = 47.867 * AMUtoKG
	TritiumMass    = 3.016 * AMUtoKG
	HydrogenMass   = 1.008 * AMUtoKG
	electronMass   = 0.00054858 * AMUtoKG // Electron mass in amu
	I              = 16.0 * eVtoJ         // Mean excitation potential in eV
	speedOfLight   = 2.997e8              // Speed of light in m/s
	eps0           = 8.8541878e-12        // Permittivity constant F/m
	electronCharge = 1.602e-19            // Coulombs
	kCoulomb       = 8.9875517873681764e9 // Coulomb's constant in N·m²/C²
	pi             = 3.141592653589793
	kc             = 2.30709e-28  // ec^2/(4 pi e0) in Joule Meter
	threshold      = 1000 * eVtoJ // cutoff in J
	mbtom2         = 1e-31        // convert millibarn to m^2
)

// Constants for the Bosch-Hale formulas
var (
	A = [][]float64{
		{69.27, 7.454e5, 2.05e3, 52.002, 0.0},
		{5750.1, 2.5226, 0.0455, 0.0, 0.0},
		{55.576, 0.21054, -3.2638e-5, 1.4987e-9, 1.8181e-13},
		{53.701, 0.33027, -1.2706e-4, 2.9327e-8, -2.515e-12},
	}
	B = [][]float64{
		{63.8, -0.995, 6.981e-5, 1.728e-4, 0.0},
		{-3.1995e-3, -8.533e-6, 5.9014e-8, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0},
	}
	BG       = []float64{34.3827, 68.7508, 31.397, 31.397}
	Q_values = []float64{17571000, 18354000, 4032000, 3268000}
)

type Atom struct {
	name    string
	density float64
	Z       int
	A       int
	m       float64
}

type Material struct {
	atoms []Atom
	mutex sync.Mutex
}

type Particle struct {
	position         float64
	energy           float64
	Z                int
	A                int
	m                float64
	scatteringEvents int
	fusionReaction   int
}

func initializeMaterial() *Material {
	return &Material{
		atoms: []Atom{
			{name: "deuterium", density: 0, Z: 1, A: 2, m: 2.014 * AMUtoKG},
			{name: "titanium", density: 5.67e28, Z: 22, A: 48, m: 47.867 * AMUtoKG},
			{name: "tritium", density: 5.67e28, Z: 1, A: 3, m: 3.016 * AMUtoKG},
			{name: "hydrogen", density: 0.0, Z: 1, A: 1, m: 1.008 * AMUtoKG},
			{name: "helium3", density: 0.0, Z: 2, A: 3, m: 3.016 * AMUtoKG},
		},
	}
}

func initializeParticles(n int) []*Particle {
	particles := make([]*Particle, n)
	for i := range particles {
		particles[i] = &Particle{position: 0.0, energy: 150000.0 * eVtoJ, Z: 1, A: 2, m: DeuteriumMass, scatteringEvents: 0, fusionReaction: -1}
	}
	return particles
}

func totalDensity(material *Material) float64 {
	total := 0.0
	for _, atom := range material.atoms {
		total += atom.density
	}
	return total
}

// Choose a particle from the material based on number density weighted probability
func chooseInteractionAtom(totalDensity float64, material *Material) Atom {
	// Calculate cumulative probabilities
	cumulativeProbabilities := make([]float64, len(material.atoms))
	cumulative := 0.0
	for i, atom := range material.atoms {
		cumulative += atom.density / totalDensity
		cumulativeProbabilities[i] = cumulative
	}

	// Generate a random number between 0 and 1
	r := rand.Float64()

	// Determine which particle to choose based on the random number
	for i, cumulativeProbability := range cumulativeProbabilities {
		if r < cumulativeProbability {
			return material.atoms[i]
		}
	}

	// Fallback (shouldn't happen if probabilities are correctly normalized)
	return material.atoms[len(material.atoms)-1]
}

func computeLosses(particle *Particle, material *Material) float64 {
	avgLoss := 0.0

	for _, atom := range material.atoms {
		//deBethe := 0.0
		deLindhardt := 0.0
		//eIz := 0.0

		// Compute Bethe energy losses
		// ionization potential eIz computed in Joules
		//if atom.Z < 13 {
		//	eIz = (12 + (7 / float64(atom.Z))) * float64(atom.Z) * eVtoJ
		//} else {
		//	eIz = (9.76 + (58.5 * math.Pow(float64(atom.Z), -1.19))) * float64(atom.Z) * eVtoJ
		//}
		// From wikipedia bethe bloch (relativistic)
		//BetaSq := (math.Pow(particle.energy, 2) + 2*math.Pow(speedOfLight, 2)*particle.energy*particle.m) / math.Pow(particle.energy+math.Pow(speedOfLight, 2)*particle.m, 2)
		//newBethe := -4 * pi / (electronMass * math.Pow(speedOfLight, 2)) * atom.density * float64(atom.Z) * math.Pow(float64(particle.Z), 2) / BetaSq * math.Pow(kc, 2) * (math.Log(2*electronMass*math.Pow(speedOfLight, 2)*BetaSq/(eIz*(1-BetaSq)) - BetaSq))

		// Bethe bloch From Nastasi 5.38.  is only positive for Deuterium incident on metals at more than ~200 keV (for titanium) and ~600 for Pd
		// logVal := math.Log(4 * particle.energy * electronMass / (particle.m * eIz))
		// deBethe = 2 * pi * float64(particle.Z*particle.Z) * float64(atom.Z) * atom.density * math.Pow(kc, 2) * (particle.m / electronMass) * logVal / particle.energy

		// Compute Lindhardt energy losses (comes out in eV/cm, must be converted to J/m)
		// k comes out dimensionless, but M1 must be in AMU and E1 in keV. Also, the density is converted to cm^3
		Se := 3.83 * math.Pow(float64(particle.Z), 7.0/6.0) * float64(atom.Z) / math.Pow(math.Pow(float64(particle.Z), 2.0/3.0)+math.Pow(float64(atom.Z), 2.0/3.0), 1.5) * math.Sqrt(particle.energy*(JtoeV/1000)/(particle.m/AMUtoKG))
		deLindhardt = Se * 1e-15 * atom.density / 1e6 * eVtoJ * 100 //1e-15 is from nastasi, 1e6 is density to cm^3, and eVtoJ*100 is to convert from eV/cm to J/m

		// combine the energy losses for a smooth transition low to high energy
		// var deComp float64
		// if deBethe > 0 {
		//	deComp = deLindhardt * deBethe / (deLindhardt + deBethe)
		//} else {
		//	deComp = deLindhardt
		//}
		avgLoss += atom.density / totalDensity(material) * deLindhardt
	}

	return avgLoss
}

func chooseReaction(particle *Particle, interactionAtom *Atom) int {
	// Check if both projectile and target are deuterium
	if particle.Z == 1 && particle.A == 2 && interactionAtom.Z == 1 && interactionAtom.A == 2 {
		// This is a DD reaction. Select between T (reaction 2) and He-3 (reaction 3) branches randomly with equal probability
		if rand.Intn(2) == 0 {
			return 2
		}
		return 3
	}
	// Projectile is deuterium. Check if target is He-3
	if particle.Z == 1 && particle.A == 2 {
		if interactionAtom.Z == 2 && interactionAtom.A == 3 {
			// If it is, select reaction 1
			return 1
		}
		// Check if target is tritium
		if interactionAtom.Z == 1 && interactionAtom.A == 3 {
			// If it is, select reaction 0
			return 0
		}
		// If target is other, return 0
		return 0
	}
	// Target is deuterium. Check if projectile is He-3
	if interactionAtom.Z == 1 && interactionAtom.A == 2 {
		if particle.Z == 2 && particle.A == 3 {
			// If it is, select reaction 1
			return 1
		}
		// Check if projectile is tritium
		if particle.Z == 1 && particle.A == 3 {
			// If it is, select reaction 0
			return 0
		}
		// If projectile is not tritium or He-3, no fusion will happen
		return 0
	}
	// If neither target nor projectile are deuterium, there are no fusion reaction
	return 0
}

func fusionCrossSection(particle *Particle, interactionAtom *Atom) float64 {
	currentReaction := chooseReaction(particle, interactionAtom)

	eCMKev := interactionAtom.m / (interactionAtom.m + particle.m) * particle.energy * JtoeV / 1000
	dummy1 := 0.0
	dummy2 := 1.0

	// Calculate the astrophysical factor S
	for i := 0; i < 4; i++ {
		dummy1 += A[currentReaction][i] * math.Pow(eCMKev, float64(i))
		dummy2 += B[currentReaction][i] * math.Pow(eCMKev, float64(i+1))
	}
	sFus := dummy1 / dummy2
	particle.fusionReaction = currentReaction

	//fmt.Printf("fusionCS: %.5e", sFus*1e3/math.Exp(BG[currentReaction]/math.Sqrt(eKev)))
	// Calculate and return the fusion cross-section in millibarn
	return sFus * 1e3 / (eCMKev * math.Exp(BG[currentReaction]/math.Sqrt(eCMKev))) * mbtom2
}

func scatteringCrossSection(interactionAtom *Atom, energy float64, material *Material) float64 {
	return 2.5e-22
}

func meanFreePath(scatteringCS float64, fusionCS float64, material *Material) float64 {
	totalDensity := totalDensity(material)

	return 1.0 / ((scatteringCS + fusionCS) * totalDensity)
}

func runSimulation(nParticles int, wg *sync.WaitGroup) []*Particle {
	defer wg.Done()

	particleStack := make([]*Particle, nParticles)
	material := initializeMaterial()
	particles := initializeParticles(nParticles)
	total_density := totalDensity(material)
	fusions := 0

	var particleWg sync.WaitGroup
	particleWg.Add(nParticles)

	for i := 0; i < nParticles; i++ {
		go func(particle *Particle, id int) {
			defer particleWg.Done()
			//fmt.Printf("Goroutine %d started, particle energy %f\n", id, particle.energy) // Print the start of the goroutine
			for particle.energy > threshold {

				// These don't need updated material parameters
				interactionAtom := chooseInteractionAtom(total_density, material)
				scatteringCS := scatteringCrossSection(&interactionAtom, particle.energy, material)

				// fusion needs updated material parameters
				material.mutex.Lock()
				fusionCS := 0.0
				if interactionAtom.Z < 3 {
					fusionCS = fusionCrossSection(particle, &interactionAtom)
				}

				//fmt.Printf("energy: %f, fusion CS: %0.5e\n", particle.energy*JtoeV/1000, fusionCS/mbtom2)
				// check if fusion or scattering occurs
				if rand.Float64() < fusionCS/(fusionCS+scatteringCS) {
					// fusion happened
					for i := range material.atoms {
						switch material.atoms[i].name {
						case "tritium":
							material.atoms[i].density += 0.1
						case "hydrogen":
							material.atoms[i].density += 0.1
						case "deuterium":
							material.atoms[i].density -= 0.1
						}
					}
					//fmt.Printf("fusion happened\n")
					fusions += 1
					particle.energy = 0.0
				} else {
					// scattering happened
					//fmt.Printf("scattering happened\n")
					//reset fusionreaction to -1
					particle.fusionReaction = -1
					particle.scatteringEvents += 1
				}
				material.mutex.Unlock()

				//fmt.Printf("scatCS %.5e, fusCS %.5e\n", scatteringCS, fusionCS)
				// Compute new location and energy
				meanPath := meanFreePath(scatteringCS, fusionCS, material)
				// compute the step length change
				stepLength := -meanPath * math.Log(rand.Float64())
				losses := computeLosses(particle, material)
				particle.energy -= losses * stepLength
				particle.position += stepLength

				// Print the particle's state
				//fmt.Printf("Goroutine %d, Particle energy(keV): %f, position(nm): %0.3e, losses(keV/nm): %f, steplength(nm): %f\n", id, particle.energy*JtoeV/1000, particle.position*1e9, losses*JtoeV/1000*1e-9, stepLength*1e9)
			}

			particleStack[id] = particle

		}(particles[i], i)
	}

	particleWg.Wait()
	fmt.Printf("fusion events %d\n", fusions)

	return particleStack
}

func main() {
	start := time.Now()
	//rand.Seed(time.Now().UnixNano())
	var wg sync.WaitGroup
	nParticles := 2_000_000

	wg.Add(1)
	particleStack := runSimulation(nParticles, &wg) // Run the simulation and get the particle stack

	wg.Wait()
	duration := time.Since(start)
	fmt.Printf("Total simulation time: %s\n", duration)

	// Write the particle stack to a CSV file
	filename := "particles.csv"
	if err := writeParticlesToCSV(particleStack, filename); err != nil {
		fmt.Printf("Error writing to CSV: %v\n", err)
		return
	}

	fmt.Printf("Particle data written to %s\n", filename)
}

func writeParticlesToCSV(particles []*Particle, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write CSV header
	header := []string{"Particle ID", "Energy (keV)", "Position (nm)", "Scattering Events", "Fusion Reactions"}
	if err := writer.Write(header); err != nil {
		return err
	}

	// Write particle data
	for i, particle := range particles {
		record := []string{
			strconv.Itoa(i),
			fmt.Sprintf("%0.2e", particle.energy*JtoeV/1000),
			fmt.Sprintf("%0.2e", particle.position*1e9),
			strconv.Itoa(particle.scatteringEvents),
			strconv.Itoa(particle.fusionReaction),
		}
		if err := writer.Write(record); err != nil {
			return err
		}
	}

	return nil
}
