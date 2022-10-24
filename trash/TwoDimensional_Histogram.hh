#ifndef TWODIMENSIONAL_HISTOGRAM_HH
#define TWODIMENSIONAL_HISTOGRAM_HH

#include "globals.hh"
#include <vector>
#include "G4SystemOfUnits.hh"

typedef std::vector<G4double> G4double_vector;



class TwoDimensional_Histogram
{
private:
    TwoDimensional_Histogram(int N_cols,      G4double col_min_value , G4double col_max_value,
                              int N_rows,      G4double row_min_value , G4double row_max_value,
                              int N_data_bins, G4double data_min_value, G4double data_max_value);
public:
   // TwoDimensional_Histogram();

    static TwoDimensional_Histogram* instance;

    int N_cols;
    G4double col_min_value;
    G4double col_max_value;
    G4double col_value_width;

    int N_rows;
    G4double row_min_value;
    G4double row_max_value;
    G4double row_value_width;

    int N_data_bins;
    G4double data_min_value;
    G4double data_max_value;
    G4double data_bin_width;

    int N_total;
    G4double_vector histogram_matrix;
    G4double_vector columns_map;
    G4double_vector rows_map;
    G4double_vector data_map;


    ~TwoDimensional_Histogram() {};

    inline TwoDimensional_Histogram* getInstance(){return instance;};

//    void init_size(int N_cols, G4double col_min_value, G4double col_max_value,
//                    int N_rows, G4double row_min_value, G4double row_max_value,
//                    int N_data_bins, G4double data_min_value,G4double data_max_value);



    void accumulate(G4double row_value, G4double col_value, G4double data_value);

    void accumulate(int row_bin_Num, G4double col_bin_Num, G4double data_bin_Num);

private:
    bool isInRange(G4double value, G4double min_value, G4double max_value, G4double value_step);
    int getBinNumber(G4double value, G4double min_value, G4double max_value, G4double value_step);
    inline int getIndex(int col_bin_num,int row_bin_num,int data_bin_num){
        return col_bin_num + (row_bin_num  + data_bin_num * N_rows) * N_cols;
    };

};

#ifdef PhaseSpaceClass_HH
#define PhaseSpaceClass_HH
class PhaseSpace{
private:
    // PhaseSpace is a singleton object (1 instance only)
    PhaseSpace(){};
public:
    static PhaseSpace* instance;
    unsigned int Nx, Ny;
    TwoDimensional_Histogram *energy; 	 // pointer to Energy matrix
    TwoDimensional_Histogram *polar_angle; 	 // pointer to polar angle matrix
    TwoDimensional_Histogram *azimuthal_angle; 	 // pointer to azimuthal angle matrix
    TwoDimensional_Histogram *fluence;   // pointer to fluence matrix

    ~PhaseSpace(){};
    static PhaseSpace* GetInstance(){
        if (instance == 0) instance = new PhaseSpace;
        return instance;
    };

};
#endif

#endif // TWODIMENSIONAL_HISTOGRAM_HH
