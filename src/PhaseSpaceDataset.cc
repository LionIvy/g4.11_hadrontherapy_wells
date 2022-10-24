#include "PhaseSpaceDataset.hh"


//======================================================================================================================================
//                          PhaseSpaceDataset
//======================================================================================================================================

PhaseSpaceDataset* PhaseSpaceDataset::instance = nullptr;
PhaseSpaceDataset::PhaseSpaceDataset(int new_N_cols, G4double new_col_min_value, G4double new_col_max_value, // X coordinate
                                     int new_N_rows, G4double new_row_min_value, G4double new_row_max_value, // Y coordinate
                                     int new_N_energy_bins, G4double new_energy_min_value, G4double new_energy_max_value) // Energy dataset
{
        G4double temp = 0;
        this->N_cols = new_N_cols;
        this->col_min_value = new_col_min_value;
        this->col_max_value = new_col_max_value;
        this->col_value_width = (col_max_value-col_min_value)/N_cols;
        //G4double_vector columns_map;
        temp = col_min_value;
        columns_map.clear();
        while(temp < col_max_value){
            columns_map.push_back(temp);
            temp += col_value_width;
        }

        this->N_rows = new_N_rows;
        this->row_min_value = new_row_min_value;
        this->row_max_value = new_row_max_value;
        this->row_value_width = (row_max_value-row_min_value)/N_rows;
        //G4double_vector rows_map;
        temp = row_min_value;
        rows_map.clear();
        while(temp < row_max_value){
            rows_map.push_back(temp);
            temp += row_value_width;
        }


        //Область флюенса
        N_fluence_total = N_rows * N_cols;

        // Область энергий
        this->N_energy_bins = new_N_energy_bins;
        this->energy_min_value = new_energy_min_value;
        this->energy_max_value = new_energy_max_value;
        this->energy_bin_width = (energy_max_value-energy_min_value)/N_energy_bins;
        //G4double_vector energy_value_map;
        energy_value_map.clear();
        temp =energy_min_value;
        while(temp < energy_max_value){
            energy_value_map.push_back(temp);
            temp += energy_bin_width;
        }
        N_energy_total = N_fluence_total * N_energy_bins;
        G4cout << "PhaseSpaceDataset constructed" << G4endl;

}
PhaseSpaceDataset::~PhaseSpaceDataset(){

}
PhaseSpaceDataset* PhaseSpaceDataset::getInstance(){
    return instance;
}
PhaseSpaceDataset* PhaseSpaceDataset::getInstance(int N_cols, G4double col_min_value, G4double col_max_value, // X coordinate
                                                  int N_rows, G4double row_min_value, G4double row_max_value, // Y coordinate
                                                  int N_energy_bins, G4double energy_min_value, G4double energy_max_value)
{

    G4cout << "if (instance) delete instance;" << G4endl;
    if (instance) delete instance;
    G4cout << "instance = new PhaseSpaceDataset" << G4endl;
    instance = new PhaseSpaceDataset(N_cols, col_min_value, col_max_value, // X coordinate
                                     N_rows, row_min_value, row_max_value, // Y coordinate
                                     N_energy_bins, energy_min_value, energy_max_value);
    instance ->initialize();
    return instance;
}


//{

//}

void PhaseSpaceDataset::initialize(){
    G4cout << "PhaseSpaceDataset initialize" << G4endl;
    energy = new TwoDimensionalHistogram(N_cols, col_min_value, col_max_value,
                               N_rows, row_min_value, row_max_value,
                               N_energy_bins, energy_min_value, energy_max_value);

    fluence = new TwoDimensionalHistogram(N_cols, col_min_value, col_max_value,
                               N_rows, row_min_value, row_max_value,
                               1, 0, 1);
    G4cout << "PhaseSpaceDataset initialized" << G4endl;
}

//Fill energy dataset
void PhaseSpaceDataset::fill_energy_dataset(G4double X, G4double Y, G4double energy_value)
{
    energy->accumulate(X, Y, energy_value);
}
void PhaseSpaceDataset::fill_energy_dataset(int Detector_Copy_Num, G4double energy_value){
    energy->accumulate(Detector_Copy_Num, energy_value);
}

//Fill fluence dataset
void PhaseSpaceDataset::fill_fluence_dataset(G4double X, G4double Y){
    fluence->accumulate(X, Y, 0.5);
}
void PhaseSpaceDataset::fill_fluence_dataset(int Detector_Copy_Num){
    fluence->accumulate(Detector_Copy_Num, 0.5);
}
