package main

import (
	"fmt"
	"math"
	"math/rand"
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
	kc             = 2.30709e-28 // ec^2/(4 pi e0) in Joule Meter
	threshold      = 500 * eVtoJ // cutoff in J
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
	position          float64
	energy            float64
	Z                 int
	m                 float64
	scattering_events int
}

func initializeMaterial() *Material {
	return &Material{
		atoms: []Atom{
			{name: "deuterium", density: 5.67e28, Z: 1, A: 2, m: 2.014 * AMUtoKG},
			{name: "titanium", density: 5.67e28, Z: 22, A: 48, m: 47.867 * AMUtoKG},
			{name: "tritium", density: 0.0, Z: 1, A: 3, m: 3.016 * AMUtoKG},
			{name: "hydrogen", density: 0.0, Z: 1, A: 1, m: 1.008 * AMUtoKG},
			{name: "helium3", density: 0.0, Z: 2, A: 3, m: 3.016 * AMUtoKG},
		},
	}
}

func initializeParticles(n int) []*Particle {
	particles := make([]*Particle, n)
	for i := range particles {
		particles[i] = &Particle{position: 0.0, energy: 2000.0 * eVtoJ, Z: 1, m: DeuteriumMass, scattering_events: 0}
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
		deBethe := 0.0
		deLindhardt := 0.0
		eIz := 0.0

		// Compute Bethe energy losses
		// ionization potential eIz computed in Joules
		if atom.Z < 13 {
			eIz = (12 + (7 / float64(atom.Z))) * float64(atom.Z) * eVtoJ
		} else {
			eIz = (9.76 + (58.5 * math.Pow(float64(atom.Z), -1.19))) * float64(atom.Z) * eVtoJ
		}
		logVal := math.Log(4 * particle.energy * electronMass / (particle.m * eIz))
		deBethe = -2 * pi * float64(particle.Z*particle.Z) * float64(atom.Z) * atom.density * math.Pow(kc, 2) * (particle.m / electronMass) * logVal / particle.energy

		fmt.Printf("energy: %f, dEdx bethe: %f\n", particle.energy*JtoeV, deBethe*JtoeV*1e-9)

		// Compute Lindhardt energy losses
		k := 3.83 * math.Pow(float64(particle.Z), 7.0/6.0) * float64(atom.Z) * atom.density / (math.Pow(particle.m, 0.5) * math.Pow(math.Pow(float64(particle.Z), 2.0/3.0)+math.Pow(float64(atom.Z), 2.0/3.0), 1.5))
		deLindhardt = 1e-23 * k * math.Sqrt(particle.energy/1000)

		fmt.Printf("energy: %f, dEdx lin: %f\n", particle.energy/eVtoJ, deLindhardt/eVtoJ*1e-9)

		// Print all values in scientific notation
		fmt.Printf("particle.energy: %.5e eV\n", particle.energy/eVtoJ)
		fmt.Printf("electronMass: %.5e kg\n", electronMass)
		fmt.Printf("deuteron mass: %.5e kg\n", particle.m)
		fmt.Printf("eIz: %.5e eV\n", eIz*JtoeV)
		fmt.Printf("logVal: %.5e\n", logVal)
		fmt.Printf("particle.Z: %d\n", particle.Z)
		fmt.Printf("atom.Z: %d\n", atom.Z)
		fmt.Printf("atom.density: %.5e atoms/m³\n", atom.density)
		fmt.Printf("kCoulomb: %.5e N·m²/C²\n", kCoulomb)
		fmt.Printf("deBethe: %.5e eV/nm\n\n", deBethe/eVtoJ*1e-9)

		var deComp float64
		if deBethe > 0 {
			deComp = deLindhardt * deBethe / (deLindhardt + deBethe)
		} else {
			deComp = deLindhardt
		}
		avgLoss += atom.density / totalDensity(material) * deComp
	}

	return avgLoss
}

func fusionCrossSection(interactionAtom *Atom, energy float64, material *Material) float64 {
	return 1e-26 * energy
}

func scatteringCrossSection(interactionAtom *Atom, energy float64, material *Material) float64 {
	return 1e-23 * energy
}

func meanFreePath(scatteringCS float64, fusionCS float64, particle *Particle, material *Material) float64 {
	totalDensity := totalDensity(material)

	return 1.0 / ((scatteringCS + fusionCS) * totalDensity)
}

func runSimulation(nParticles int, wg *sync.WaitGroup) {
	defer wg.Done()

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
				fusionCS := fusionCrossSection(&interactionAtom, particle.energy, material)

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
					particle.scattering_events += 1
				}
				material.mutex.Unlock()

				// Compute new location and energy
				meanPath := meanFreePath(scatteringCS, fusionCS, particle, material)
				losses := computeLosses(particle, material)
				particle.energy -= losses
				particle.position += meanPath

				// Print the particle's state
				//fmt.Printf("Goroutine %d, Particle energy: %f, position: %f, ses: %d\n", id, particle.energy, particle.position, particle.scattering_events)
			}

		}(particles[i], i)
	}

	particleWg.Wait()
	fmt.Printf("fusion events %d\n", fusions)
}

func main() {
	start := time.Now()
	rand.Seed(time.Now().UnixNano())
	var wg sync.WaitGroup
	nParticles := 1

	wg.Add(1)
	go runSimulation(nParticles, &wg)

	wg.Wait()
	duration := time.Since(start)
	fmt.Printf("Total simulation time: %s\n", duration)
}
