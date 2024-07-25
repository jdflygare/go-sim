package main

import (
	"fmt"
	"math/rand"
	"time"
)

// Choose a number based on given probabilities
func choose(numbers []float64, probabilities []float64) float64 {
	if len(numbers) != len(probabilities) {
		panic("Numbers and probabilities must have the same length")
	}

	// Calculate cumulative probabilities
	cumulativeProbabilities := make([]float64, len(probabilities))
	cumulative := 0.0
	for i, probability := range probabilities {
		cumulative += probability
		cumulativeProbabilities[i] = cumulative
	}

	// Generate a random number between 0 and 1
	r := rand.Float64()

	// Determine which number to choose based on the random number
	for i, cumulativeProbability := range cumulativeProbabilities {
		if r < cumulativeProbability {
			return numbers[i]
		}
	}

	// Fallback (shouldn't happen if probabilities are correctly normalized)
	return numbers[len(numbers)-1]
}

func main() {
	rand.Seed(time.Now().UnixNano())

	numbers := []float64{1, 2, 3}
	probabilities := []float64{1.0, 0.01, 0.03} // Example probabilities

	chosen := choose(numbers, probabilities)

	fmt.Printf("Chosen number: %.1f\n", chosen)
}
