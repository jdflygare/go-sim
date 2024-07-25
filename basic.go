package main

import (
	"fmt"
	"math/rand"
	"sync"
	"time"
)

// All densities are number densities
// All energies are in eV

// Global constants
const (
	AvogadroNumber = 6.02214076e23     // Avogadro's number
	AMUtoG         = 1.66053906660e-24 // Conversion factor from AMU to g
	DeuteriumMass  = 2.014 * AMUtoG    // molecular weights are in g/mol
	TitaniumMass   = 47.867 * AMUtoG
	TritiumMass    = 3.016 * AMUtoG
	HydrogenMass   = 1.008 * AMUtoG
	MeVtoEV        = 1.0e6           // Conversion from MeV to eV
	electronMass   = 0.511 * MeVtoEV // Electron mass in eV
	K              = 0.307           // MeV cm^2/g
	I              = 16.0 * 1.0e-6   // Mean excitation potential in eV (typical value for materials, to be adjusted based on material)
	speedOfLight   = 3.0e8           // Speed of light in m/s
)

type Atom struct {
	name    string
	density float64
	Z       int
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
			{name: "deuterium", density: 5.67e22, Z: 1, m: 2.014 * AMUtoG},
			{name: "titanium", density: 5.67e22, Z: 22, m: 47.867 * AMUtoG},
			{name: "tritium", density: 0.0, Z: 1, m: 3.016 * AMUtoG},
			{name: "hydrogen", density: 0.0, Z: 1, m: 1.008 * AMUtoG},
			{name: "helium3", density: 0.0, Z: 2, m: 3.016 * AMUtoG},
		},
	}
}

func initializeParticles(n int) []*Particle {
	particles := make([]*Particle, n)
	for i := range particles {
		particles[i] = &Particle{position: 0.0, energy: 2000.0, Z: 1, m: DeuteriumMass, scattering_events: 0}
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
			for particle.energy > 0.01 {

				// These don't need updated material parameters
				interactionAtom := chooseInteractionAtom(total_density, material)
				scatteringCS := scatteringCrossSection(&interactionAtom, particle.energy, material)

				// fusion needs updated material parameters
				material.mutex.Lock()
				fusionCS := fusionCrossSection(particle.energy, material)

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
				meanPath := meanFreePath(particle, material)
				losses := computeLosses(particle, material)
				particle.energy -= losses
				particle.position += meanPath

				// Print the particle's state
				fmt.Printf("Goroutine %d, Particle energy: %f, position: %f, ses: %d\n", id, particle.energy, particle.position, particle.scattering_events)
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
	nParticles := 1 //1_000_00

	wg.Add(1)
	go runSimulation(nParticles, &wg)

	wg.Wait()
	duration := time.Since(start)
	fmt.Printf("Total simulation time: %s\n", duration)
}
