package main

import (
	"fmt"
	"log"
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
	kc             = 2.30709e-28    // ec^2/(4 pi e0) in Joule Meter
	threshold      = 1000 * eVtoJ   // cutoff in J
	mbtom2         = 1e-31          // convert millibarn to m^2
	kqq            = 1.43965        // keV pm
	kqqSI          = 2.30657e-28    // Coulomb^2 Meter/Farad
	sro            = 1.25e-15       // m
	srSI           = 0.00327579e-12 // m
	hbar           = 1.0545e-34     // J s
	alpha          = 0.00730013     // fine structure constant
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
	fusionReaction   []int
	fusionEnergy     float64
	enhancement      float64
}

type LookupTable struct {
	xData, yData []float64
}

func initializeMaterial() *Material {
	//ti = 5.67e28, pd=6.79e28
	return &Material{
		atoms: []Atom{
			{name: "deuterium", density: 0, Z: 1, A: 2, m: 2.014 * AMUtoKG},
			{name: "titanium", density: 5.67e28, Z: 22, A: 48, m: 47.867 * AMUtoKG},
			{name: "palladium", density: 0, Z: 46, A: 106, m: 106.42 * AMUtoKG},
			{name: "tritium", density: 5.67e28, Z: 1, A: 3, m: 3.016 * AMUtoKG},
			{name: "hydrogen", density: 0.0, Z: 1, A: 1, m: 1.008 * AMUtoKG},
			{name: "helium3", density: 0.0, Z: 2, A: 3, m: 3.016 * AMUtoKG},
		},
	}
}

func initializeParticles(n int) []*Particle {
	particles := make([]*Particle, n)
	for i := range particles {
		Ek := 100000.0
		Espread := Ek + (((rand.Float64() * 2) - 1) * (0.05 * Ek))
		particles[i] = &Particle{position: 0.0, energy: Espread * eVtoJ, Z: 1, A: 2, m: DeuteriumMass, scatteringEvents: 0, fusionReaction: []int{}, fusionEnergy: 0.0, enhancement: 1.0}
	}
	return particles
}

func runSimulation(nParticles int, replicas int, wg *sync.WaitGroup) ([]*Particle, []*Particle) {
	defer wg.Done()

	particleStack := make([]*Particle, nParticles)
	fusionStack := []*Particle{}
	material := initializeMaterial()
	particles := initializeParticles(nParticles)
	total_density := totalDensity(material)
	//var fusionStackMutex sync.Mutex
	// load scattering data and initialize lookup table
	scatData, err := LoadCSV("scat_results_5e-5_keV.csv")
	if err != nil {
		log.Fatal(err)
	}
	// load enhancement data
	enhData, err := LoadCSV("enhancements/all_enh.csv")
	if err != nil {
		log.Fatal(err)
	}

	fmt.Print("Initializing scattering tables\n")
	ltPd := NewLookupTable(scatData, 1, 2)
	ltTi := NewLookupTable(scatData, 3, 4)
	ltD := NewLookupTable(scatData, 5, 6)
	ltTr := NewLookupTable(scatData, 7, 8)
	//fmt.Print(lt)

	fmt.Print("Initializing enhancement tables\n")
	//enh110 := NewLookupTable(enhData, 0, 1)
	//enh125 := NewLookupTable(enhData, 0, 2)
	//enh300 := NewLookupTable(enhData, 0, 3)
	//enh600 := NewLookupTable(enhData, 0, 4)
	enh1000 := NewLookupTable(enhData, 0, 5)
	//enh3000 := NewLookupTable(enhData, 0, 6)
	//enh4400 := NewLookupTable(enhData, 0, 7)
	//enh40400 := NewLookupTable(enhData, 0, 8)

	fusions := 0

	var particleWg sync.WaitGroup
	particleWg.Add(nParticles)

	for i := 0; i < nParticles; i++ {
		go func(particle *Particle, id int) {
			defer particleWg.Done()

			//track replicas to not duplicate particles
			internal_replicas := replicas

			for particle.energy > threshold {

				// ********* Choose interaction partner and compute center of mass energy ***********
				interactionAtom := chooseInteractionAtom(total_density, material)
				eCMKev := interactionAtom.m / (interactionAtom.m + particle.m) * particle.energy * JtoeV / 1000
				//eCM := interactionAtom.m / (interactionAtom.m + particle.m) * particle.energy

				//  *********** compute scattering cross section ***********
				scatteringCS := 0.0
				//is returned in cm^2, convert to m^2
				if interactionAtom.Z == 46 {
					scatteringCS = ltPd.InterpolateValue(eCMKev) * 1e-6
				}
				if interactionAtom.Z == 22 {
					scatteringCS = ltTi.InterpolateValue(eCMKev) * 1e-6
				}
				if interactionAtom.A == 2 {
					scatteringCS = ltD.InterpolateValue(eCMKev) * 1e-6
				}
				if interactionAtom.A == 3 {
					scatteringCS = ltTr.InterpolateValue(eCMKev) * 1e-6
				}

				// **************  compute fusion cross section *************
				//material.mutex.Lock()
				fusionCS := 0.0
				reactionType := 0

				if interactionAtom.Z < 3 {
					fusionCS, reactionType = fusionCrossSection(eCMKev, particle, &interactionAtom)
					if particle.energy*JtoeV/1000 <= 100.0 {
						fue := enh1000.InterpolateValue(particle.energy * JtoeV / 1000)
						// include scattering effects
						fue = fue * particle.energy / (particle.energy + 1000*eVtoJ)
						//fue = 1
						particle.enhancement = fue
					}
					//fmt.Printf("fue: %0.5e\n", fue)

					fusionCS = particle.enhancement * fusionCS
				}

				// loop over this X times to boost statistics for fusion events
				cs_rat := fusionCS / (fusionCS + scatteringCS)

				for i := 0; i < internal_replicas; i++ {
					//fmt.Printf("energy: %f, fusion CS: %0.5e\n", particle.energy*JtoeV/1000, fusionCS/mbtom2)
					// check if fusion or scattering occurs

					if rand.Float64() < cs_rat && interactionAtom.Z < 3 {
						// fusion happened
						fusions += 1
						particle.fusionReaction = append(particle.fusionReaction, reactionType)
						particle.fusionEnergy = particle.energy
						particleCopy := *particle
						if particle != nil {
							fusionStack = append(fusionStack, &particleCopy)
						}

						//fmt.Printf("Particle energy(keV): %f, position(nm): %0.3e, int.A: %d, eCM(keV): %f, fue: %f\n", particle.energy*JtoeV/1000, particle.position*1e9, interactionAtom.A, eCMKev, particle.enhancement)
						//kill the replica
						if internal_replicas > 0 {
							internal_replicas -= 1
						} else {
							particle.energy = 0.0
						}
						//particle.energy = 0.0
					} else {
						// scattering happened
						particle.scatteringEvents += 1
					}
				}
				//material.mutex.Unlock()

				// Compute new location and energy
				meanPath := meanFreePath(scatteringCS, fusionCS, material, &interactionAtom)
				// compute the step length change
				stepLength := meanPath //* math.Log(rand.Float64()) //-meanPath * math.Log(rand.Float64()) * 2.2
				losses := computeLosses(particle, material)

				// Print the particle's state
				//fmt.Printf("Particle energy(keV): %f, atom: %d, CMenergy(keV): %f, scatCS(mb): %0.3e, fusCS(mb): %0.3e, losses(keV/nm): %f\n", particle.energy*JtoeV/1000, interactionAtom.Z, eCMKev, scatteringCS/mbtom2, fusionCS/mbtom2, losses*JtoeV/1000*1e-9)

				//fmt.Printf("Particle energy(keV): %f, position(nm): %0.3e, int.A: %d, eCM(keV): %f, fue: %f\n", particle.energy*JtoeV/1000, particle.position*1e9, interactionAtom.A, eCMKev, particle.enhancement)

				// Update particle energy and position
				particle.energy -= losses * stepLength
				particle.position += stepLength

			}

			particleStack[id] = particle

		}(particles[i], i)
	}

	particleWg.Wait()
	fmt.Printf("fusion events %d\n", fusions)

	return particleStack, fusionStack
}

func main() {

	start := time.Now()
	//rand.Seed(time.Now().UnixNano())
	var wg sync.WaitGroup
	nParticles := 1_000_000
	replicas := 10000

	wg.Add(1)
	particleStack, fusionStack := runSimulation(nParticles, replicas, &wg) // Run the simulation and get the particle stack

	wg.Wait()
	duration := time.Since(start)
	fmt.Printf("Total simulation time: %s,  Total particles: %0.1e, particle stack length: %d\n", duration, float64(nParticles*replicas), len(particleStack))

	// Write the particle stack to a CSV file
	filename := "/Users/josh/Documents/Fusion/py_data/outputs/TiTr_100keV_1000eV_1mx10000.csv"
	if err := writeParticlesToCSV(fusionStack, filename); err != nil {
		fmt.Printf("Error writing to CSV: %v\n", err)
		return
	}

	fmt.Printf("Particle data written to %s\n", filename)
}
