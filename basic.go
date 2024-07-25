package main

import (
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
	position float64
	energy   float64
	Z        int
	m        float64
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
		particles[i] = &Particle{position: 0.0, energy: 1.0, Z: 1, m: DeuteriumMass}
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
	for i, particle := range material.atoms {
		cumulative += particle.density / totalDensity
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

func computeLosses(energy float64, material *Material) float64 {
	return energy * 0.05
}

func fusionCrossSection(interactionParticle *Atom, energy float64, material *Material) float64 {
	return 1e-24 * energy
}

func scatteringCrossSection(interactionParticle *Atom, energy float64, material *Material) float64 {
	return 1e-25 * energy
}

func meanFreePath(sigmaf float64, sigmas float64, particle *Particle, material *Material) float64 {

	return 1.0 / ((sigmaf + sigmas) * totalDensity(material))
}

func runSimulation(nParticles int, wg *sync.WaitGroup) {
	defer wg.Done()

	material := initializeMaterial()
	particles := initializeParticles(nParticles)
	total_density := totalDensity(material)

	var particleWg sync.WaitGroup
	particleWg.Add(nParticles)

	for i := 0; i < nParticles; i++ {
		go func(particle *Particle, id int) {
			defer particleWg.Done()
			for particle.energy > 0.01 {

				// These don't need updated material parameters
				interactionParticle = chooseInteractionParticle(total_density, material)
				scatteringCS := scatteringCrossSection(ineractionParticle, particle.energy, material)
				//fmt.Printf("Goroutine: %d, Material before lock and loop: %+v\n", id, material)

				// fusion needs updated material parameters
				material.mutex.Lock()
				fusionCS := fusionCrossSection(particle.energy, material)

				if rand.Float64() < fusionCS/(fusionCS+scatteringCS) {
					//material.mutex.Lock()
					material.tritiumDensity += 0.1
					material.hydrogenDensity += 0.1
					material.deuteriumDensity -= 0.1
					//fmt.Printf("Goroutine: %d, Material after fusion: %+v\n", id, material)
					//material.mutex.Unlock()
					beam_particle.energy = 0.0
				} else {
					beam_particle.energy *= 0.9
					//material.mutex.Lock()
					material.tritiumDensity += 0.1
					material.hydrogenDensity += 0.1
					material.deuteriumDensity -= 0.1
					//fmt.Printf("Goroutine: %d, Material after scattering: %+v\n", id, material)
					//material.mutex.Unlock()
				}
				material.mutex.Unlock()
				//fmt.Printf("Goroutine: %d, Material after loop: %+v\n", id, material)

				// compute new location and energy
				meanPath := meanFreePath(particle, material)
				losses := computeLosses(particle.energy, material)
				particle.energy -= losses
				particle.position += meanPath
			}

		}(particles[i], i)
	}

	particleWg.Wait()
}

func main() {
	rand.Seed(time.Now().UnixNano())
	var wg sync.WaitGroup
	nParticles := 1_000_000

	wg.Add(1)
	go runSimulation(nParticles, &wg)

	wg.Wait()
}
