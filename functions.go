package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"
)

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
		// From wikipedia bethe bloch (relativistic)
		//BetaSq := (math.Pow(particle.energy, 2) + 2*math.Pow(speedOfLight, 2)*particle.energy*particle.m) / math.Pow(particle.energy+math.Pow(speedOfLight, 2)*particle.m, 2)
		//newBethe := -4 * pi / (electronMass * math.Pow(speedOfLight, 2)) * atom.density * float64(atom.Z) * math.Pow(float64(particle.Z), 2) / BetaSq * math.Pow(kc, 2) * (math.Log(2*electronMass*math.Pow(speedOfLight, 2)*BetaSq/(eIz*(1-BetaSq)) - BetaSq))

		// Bethe bloch From Nastasi 5.38.  is only positive for Deuterium incident on metals at more than ~200 keV (for titanium) and ~600 for Pd
		logVal := math.Log(4 * particle.energy * electronMass / (particle.m * eIz))
		deBethe = 2 * pi * float64(particle.Z*particle.Z) * float64(atom.Z) * atom.density * math.Pow(kc, 2) * (particle.m / electronMass) * logVal / particle.energy

		// Compute Lindhardt energy losses (comes out in eV/cm, must be converted to J/m)
		// k comes out dimensionless, but M1 must be in AMU and E1 in keV. Also, the density is converted to cm^3
		Se := 3.83 * math.Pow(float64(particle.Z), 7.0/6.0) * float64(atom.Z) / math.Pow(math.Pow(float64(particle.Z), 2.0/3.0)+math.Pow(float64(atom.Z), 2.0/3.0), 1.5) * math.Sqrt(particle.energy*(JtoeV/1000)/(particle.m/AMUtoKG))
		deLindhardt = Se * 1e-15 * atom.density / 1e6 * eVtoJ * 100 //1e-15 is from nastasi, 1e6 is density to cm^3, and eVtoJ*100 is to convert from eV/cm to J/m

		// combine the energy losses for a smooth transition low to high energy
		var deComp float64
		if deBethe > 0 {
			deComp = deLindhardt * deBethe / (deLindhardt + deBethe)
		} else {
			deComp = deLindhardt
		}
		avgLoss += deComp
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

func fusionCrossSection(eCMKev float64, particle *Particle, interactionAtom *Atom) (float64, int) {
	currentReaction := chooseReaction(particle, interactionAtom)

	dummy1 := 0.0
	dummy2 := 1.0

	// Calculate the astrophysical factor S
	for i := 0; i < 4; i++ {
		dummy1 += A[currentReaction][i] * math.Pow(eCMKev, float64(i))
		dummy2 += B[currentReaction][i] * math.Pow(eCMKev, float64(i+1))
	}
	sFus := dummy1 / dummy2

	//fmt.Printf("fusionCS: %.5e", sFus*1e3/math.Exp(BG[currentReaction]/math.Sqrt(eKev)))
	// Calculate and return the fusion cross-section in m^2
	sigmaFus := sFus * 1e3 / (eCMKev * math.Exp(BG[currentReaction]/math.Sqrt(eCMKev))) * mbtom2

	return sigmaFus, currentReaction
}

func getScreeningEnhancement(particle *Particle, interactionAtom *Atom, ECM float64, Ue float64) float64 {
	// Reduced mass (mu)
	mu := interactionAtom.m * particle.m / (interactionAtom.m + particle.m)

	// Strong radius position (sr)
	sr := sro * (math.Pow(float64(particle.A), 1.0/3.0) + math.Pow(float64(interactionAtom.A), 1.0/3.0))

	// Coulomb potential (Vc)
	Vc := kqqSI * float64(interactionAtom.Z) / sr

	// Gamow energy (Eg)
	Eg := 2 * mu * math.Pow(speedOfLight, 2.0) * math.Pow(math.Pi*alpha*float64(particle.Z)*float64(interactionAtom.Z), 2.0)

	// Compute the Gamow factors
	Gcbase := GamowC(Eg, ECM, Vc)
	//fmt.Printf("Gcbase %f\n", Gcbase)

	Gcenh := GamowC(Eg, ECM+Ue, Vc)
	//fmt.Printf("GcEnh %f\n", Gcenh)

	// Compute the screening enhancement factor
	fue := ECM / (ECM + Ue) * math.Exp(Gcbase-Gcenh)

	return fue
}

func GamowC(Eg float64, ECM float64, Vc float64) float64 {
	Gc := math.Pow(Eg/ECM, 0.5) * (2 / math.Pi * 1 / math.Cos(math.Sqrt(ECM/Vc)-math.Sqrt(ECM/Vc*(1-ECM/Vc))))
	return Gc
}

func scatteringCrossSection(interactionAtom *Atom, energyCM float64) float64 {
	scatCS := 0.0
	// energy comes in keV center of mass, and function returns cm^2

	if interactionAtom.Z == 46 {
		scatCS = 2.24297e-13 + 1.27243e-15*math.Pow(energyCM, 1.0/3.0) - 2.18116e-13*math.Pow(energyCM, 0.01)
	}
	if interactionAtom.Z == 22 {
		scatCS = 1.5694e-12 + 9.02957e-13*math.Pow(energyCM, 1.0/40.0) - 2.46483e-12*math.Pow(energyCM, 0.01)
	}
	if interactionAtom.A == 2 {
		scatCS = 1.9747e-12 + 1.18019e-12*math.Pow(energyCM, 1.0/40.0) - 3.14955e-12*math.Pow(energyCM, 0.01)
	}
	if interactionAtom.A == 3 {
		scatCS = 2.75293e-12 + 3.82808e-12*math.Pow(energyCM, 1.0/60.0) - 6.57569e-12*math.Pow(energyCM, 0.01)
	}
	return scatCS * 1e-6
}

func meanFreePath(scatteringCS float64, fusionCS float64, material *Material, interactionAtom *Atom) float64 {
	//totalDensity := totalDensity(material)
	// scattering is always higher than fusion, so determines the shorter path.
	return 1.0 / (scatteringCS * interactionAtom.density)
}

func writeParticlesToCSV(particles []*Particle, filename string) error {
	// Create the original file with all particles
	if err := saveToFile(particles, filename); err != nil {
		return err
	}

	// Filter particles to only include those with fusionEnergy > 0
	var fusion_particles []*Particle
	for _, particle := range particles {
		if particle.fusionEnergy > 0 {
			fusion_particles = append(fusion_particles, particle)
		}
	}

	// Modify the filename to append "_fusion" before the extension for the filtered file
	if strings.HasSuffix(filename, ".csv") {
		filename = strings.TrimSuffix(filename, ".csv") + "_fusion.csv"
	} else {
		filename += "_fusion"
	}

	// Create the filtered file with only fusion_particles
	if err := saveToFile(fusion_particles, filename); err != nil {
		return err
	}

	return nil
}

func saveToFile(particles []*Particle, filename string) error {
	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write CSV header
	header := []string{"Particle ID", "Energy (keV)", "Position (nm)", "Scattering Events", "Fusion Reactions", "Fusion Energy", "Enhancement"}
	if err := writer.Write(header); err != nil {
		return err
	}

	// Write particle data
	for i, particle := range particles {
		// Convert fusionReaction slice to a single string
		fusionReactionStr := "1"
		if len(particle.fusionReaction) == 0 {
			fusionReactionStr = "-1"
		} else {
			for _, num := range particle.fusionReaction {
				fusionReactionStr += strconv.Itoa(num)
			}
		}

		record := []string{
			strconv.Itoa(i),
			fmt.Sprintf("%f", particle.energy*JtoeV/1000),
			fmt.Sprintf("%f", particle.position*1e9),
			strconv.Itoa(particle.scatteringEvents),
			fusionReactionStr, // Use the concatenated string
			fmt.Sprintf("%f", particle.fusionEnergy*JtoeV/1000),
			fmt.Sprintf("%f", particle.enhancement),
		}
		if err := writer.Write(record); err != nil {
			return err
		}
	}

	return nil
}

// LoadCSV loads the CSV data from a file.
func LoadCSV(filename string) ([][]float64, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		return nil, err
	}

	// Convert records to a slice of float64 slices
	var data [][]float64
	for _, record := range records[1:] { // skip header
		var row []float64
		for _, value := range record {
			floatValue, err := strconv.ParseFloat(value, 64)
			if err != nil {
				return nil, err
			}
			row = append(row, floatValue)
		}
		data = append(data, row)
	}
	return data, nil
}

// Interpolate performs linear interpolation.
func Interpolate(x, x0, x1, y0, y1 float64) float64 {
	return (y0 + (y1-y0)*(x-x0)/(x1-x0))
}

// NewLookupTable creates a new LookupTable from data with specified columns for x and y.
func NewLookupTable(data [][]float64, xColumn, yColumn int) *LookupTable {
	xData := make([]float64, len(data))
	yData := make([]float64, len(data))
	for i, row := range data {
		xData[i] = row[xColumn]
		yData[i] = row[yColumn]
	}

	// Ensure the data is sorted by xData
	sort.SliceStable(data, func(i, j int) bool {
		return xData[i] < xData[j]
	})

	return &LookupTable{xData: xData, yData: yData}
}

// InterpolateValue performs the interpolation for a given x value.
func (lt *LookupTable) InterpolateValue(x float64) float64 {
	i := sort.Search(len(lt.xData), func(i int) bool { return lt.xData[i] >= x })
	//fmt.Printf("found i: %d\n", i) // Adding a newline for better logging

	// Edge cases
	if i == 0 {
		return lt.yData[0]
	}
	if i == len(lt.xData) {
		return lt.yData[len(lt.yData)-1]
	}

	return Interpolate(x, lt.xData[i-1], lt.xData[i], lt.yData[i-1], lt.yData[i])
}
