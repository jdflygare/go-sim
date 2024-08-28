package main

import (
	"fmt"
	"log"
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
		particles[i] = &Particle{position: 0.0, energy: 210000.0 * eVtoJ, Z: 1, A: 2, m: DeuteriumMass, scatteringEvents: 0, fusionReaction: []int{}, fusionEnergy: 0.0, enhancement: 0.0}
	}
	return particles
}

func runSimulation(nParticles int, replicas int, wg *sync.WaitGroup) []*Particle {
	defer wg.Done()

	particleStack := make([]*Particle, nParticles)
	material := initializeMaterial()
	particles := initializeParticles(nParticles)
	total_density := totalDensity(material)
	// load scattering data and initialize lookup table
	scatData, err := LoadCSV("scat_results_1e-1.csv")
	if err != nil {
		log.Fatal(err)
	}

	fmt.Print("Initializing lookup table\n")
	ltPd := NewLookupTable(scatData, 1, 2)
	ltTi := NewLookupTable(scatData, 3, 4)
	ltD := NewLookupTable(scatData, 5, 6)
	ltTr := NewLookupTable(scatData, 7, 8)
	//fmt.Print(lt)

	fusions := 0

	var particleWg sync.WaitGroup
	particleWg.Add(nParticles)

	for i := 0; i < nParticles; i++ {
		go func(particle *Particle, id int) {
			defer particleWg.Done()
			//fmt.Printf("Goroutine %d started, particle energy %0.5e\n", id, particle.energy*JtoeV/1000) // Print the start of the goroutine
			for particle.energy > threshold {

				// Choose interaction partner and compute center of mass energy
				interactionAtom := chooseInteractionAtom(total_density, material)
				eCMKev := interactionAtom.m / (interactionAtom.m + particle.m) * particle.energy * JtoeV / 1000
				eCM := interactionAtom.m / (interactionAtom.m + particle.m) * particle.energy
				//scatteringCS := scatteringCrossSection(&interactionAtom, particle.energy, material)

				//  *********** compute scattering cross section ***********
				scatteringCS := 0.0
				//is returned in cm^2, convert to m^2
				// Goes to previously interpolated functions if particle energy (lab) is < 80 keV
				if (particle.energy * JtoeV / 1000) <= 0 {
					fmt.Print("used Mathematica fit")
					scatteringCS = scatteringCrossSection(&interactionAtom, eCMKev)
				} else {
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
				}
				//fmt.Printf("scatCS: %0.5e", scatteringCS)

				// **************  compute fusion cross section *************
				//material.mutex.Lock()
				fusionCS := 0.0
				reactionType := 0
				fue := 0.0
				if interactionAtom.Z < 3 {
					fusionCS, reactionType = fusionCrossSection(eCMKev, particle, &interactionAtom)
					fue = getScreeningEnhancement(particle, &interactionAtom, eCM, 0*eVtoJ)
					particle.enhancement = fue
					//fmt.Printf("fue: %0.5e\n", fue)

					fusionCS = fue * fusionCS
				}

				// loop over this X times to boost statistics for fusion events
				for i := 0; i < replicas; i++ {
					//fmt.Printf("energy: %f, fusion CS: %0.5e\n", particle.energy*JtoeV/1000, fusionCS/mbtom2)
					// check if fusion or scattering occurs
					if rand.Float64() < fusionCS/(fusionCS+scatteringCS) {
						// fusion happened
						/*for i := range material.atoms {
							switch material.atoms[i].name {
							case "tritium":
								material.atoms[i].density += 0.1
							case "hydrogen":
								material.atoms[i].density += 0.1
							case "deuterium":
								material.atoms[i].density -= 0.1
							}
						}*/
						//fmt.Printf("fusion happened\n")
						fusions += 1
						particle.fusionReaction = append(particle.fusionReaction, reactionType)
						particle.fusionEnergy = particle.energy
						particle.energy = 0.0
					} else {
						// scattering happened
						//fmt.Printf("scattering happened\n")
						//reset fusionreaction to -1
						//particle.fusionReaction = []int{}
						particle.scatteringEvents += 1
					}
				}
				//material.mutex.Unlock()

				//fmt.Printf("scatCS %.5e, fusCS %.5e\n", scatteringCS, fusionCS)
				// Compute new location and energy
				meanPath := meanFreePath(scatteringCS, fusionCS, material)
				// compute the step length change
				stepLength := -meanPath * math.Log(rand.Float64()) * 4.5
				losses := computeLosses(particle, material)
				particle.energy -= losses * stepLength
				particle.position += stepLength

				// Print the particle's state
				//fmt.Printf("Particle energy(keV): %f, atom: %d, CMenergy(keV): %f, scatCS(mb): %0.3e, fusCS(mb): %0.3e, losses(keV/nm): %f\n", particle.energy*JtoeV/1000, interactionAtom.Z, eCMKev, scatteringCS/mbtom2, fusionCS/mbtom2, losses*JtoeV/1000*1e-9)

				//fmt.Printf("Particle energy(keV): %f, position(nm): %0.3e, losses(keV/nm): %f, steplength(nm): %f\n", particle.energy*JtoeV/1000, particle.position*1e9, losses*JtoeV/1000*1e-9, stepLength*1e9)
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
	nParticles := 1_000_000
	replicas := 10

	wg.Add(1)
	particleStack := runSimulation(nParticles, replicas, &wg) // Run the simulation and get the particle stack

	wg.Wait()
	duration := time.Since(start)
	fmt.Printf("Total simulation time: %s\n", duration)

	// Write the particle stack to a CSV file
	filename := "/Users/josh/Documents/Fusion/py_data/outputs/TiTr_test_5e-1.csv"
	if err := writeParticlesToCSV(particleStack, filename); err != nil {
		fmt.Printf("Error writing to CSV: %v\n", err)
		return
	}

	fmt.Printf("Particle data written to %s\n", filename)
}
