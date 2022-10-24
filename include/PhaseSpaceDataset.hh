#ifndef PHASESPACEDATASET_HH
#define PHASESPACEDATASET_HH

#include "globals.hh"
#include <vector>
#include "G4SystemOfUnits.hh"
#include "TwoDimensional_Histogram.hh"

class PhaseSpaceDataset
{

private:
    // PhaseSpaceDataset is supposed to be a singleton
    PhaseSpaceDataset(int new_N_cols, G4double new_col_min_value, G4double new_col_max_value, // X coordinate
                      int new_N_rows, G4double new_row_min_value, G4double new_row_max_value, // Y coordinate
                      int new_N_energy_bins, G4double new_energy_min_value, G4double new_energy_max_value); // Energy dataset
    static PhaseSpaceDataset* instance;
    std::ofstream ofs;

public:    
    // Область построения (координаты)
    int N_cols;
    G4double col_min_value;
    G4double col_max_value;
    G4double col_value_width;
    G4double_vector columns_map;

    int N_rows;
    G4double row_min_value;
    G4double row_max_value;
    G4double row_value_width;
    G4double_vector rows_map;

    //Область флюенса
    int N_fluence_total;
    TwoDimensionalHistogram *fluence;   // pointer to fluence matrix

    // Область энергий
    int N_energy_bins;
    G4double energy_min_value;
    G4double energy_max_value;
    G4double energy_bin_width;
    G4double_vector energy_value_map;
    int N_energy_total;
    TwoDimensionalHistogram *energy; 	 // pointer to Energy matrix



//    TwoDimensional_Histogram *polar_angle; 	 // pointer to polar angle matrix
//    TwoDimensional_Histogram *azimuthal_angle; 	 // pointer to azimuthal angle matrix



    ~PhaseSpaceDataset();



    // Get object instance only
    static PhaseSpaceDataset* getInstance();

    // Make & Get instance
    // Update 4 multiple datasets is needed
    static PhaseSpaceDataset* getInstance(int N_cols, G4double col_min_value, G4double col_max_value, // X coordinate
                                          int N_rows, G4double row_min_value, G4double row_max_value, // Y coordinate
                                          int N_energy_bins, G4double energy_min_value, G4double energy_max_value);



    //All elements of dataset vectors are initialize & set to zero
    void initialize();
    //void clear();

    //Fill energy dataset
    void fill_energy_dataset(G4double X, G4double Y, G4double energy);
    void fill_energy_dataset(int Detector_Copy_Num, G4double energy);

    //Fill fluence dataset
    void fill_fluence_dataset(G4double X, G4double Y);
    void fill_fluence_dataset(int Detector_Copy_Num);

};

#endif // PHASESPACEDATASET_HH
