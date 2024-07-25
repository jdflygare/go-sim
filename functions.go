package main

import "math/rand"

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
	return rand.Float64()
}

func fusionCrossSection(energy float64, material *Material) float64 {
	return 1e-26 * energy
}

func scatteringCrossSection(interactionAtom *Atom, energy float64, material *Material) float64 {
	return 1e-23 * energy
}

func meanFreePath(particle *Particle, material *Material) float64 {
	totalDensity := totalDensity(material)
	sigmaf := fusionCrossSection(particle.energy, material)
	sigmas := scatteringCrossSection(&material.atoms[0], particle.energy, material) // Assuming scattering cross section of the first atom for simplification

	return 1.0 / ((sigmaf + sigmas) * totalDensity)
}
