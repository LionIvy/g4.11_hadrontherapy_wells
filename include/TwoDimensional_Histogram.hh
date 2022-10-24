#ifndef TWODIMENSIONALHISTOGRAM_HH
#define TWODIMENSIONALHISTOGRAM_HH

#include "globals.hh"
#include <vector>
#include "G4SystemOfUnits.hh"

#include <fstream>

typedef std::vector<G4double> G4double_vector;

class TwoDimensionalHistogram
{
public:
   // TwoDimensionalHistogram();

    //TwoDimensionalHistogram* instance;



    TwoDimensionalHistogram(){};

    TwoDimensionalHistogram(int N_cols,      G4double col_min_value , G4double col_max_value,
                              int N_rows,      G4double row_min_value , G4double row_max_value,
                              int N_data_bins, G4double data_min_value, G4double data_max_value);
    ~TwoDimensionalHistogram(){};

    void set_and_initialize(int N_cols,      G4double col_min_value , G4double col_max_value,
                            int N_rows,      G4double row_min_value , G4double row_max_value,
                            int N_data_bins, G4double data_min_value, G4double data_max_value);
    void init(int N_cols,      G4double col_min_value , G4double col_max_value,
                            int N_rows,      G4double row_min_value , G4double row_max_value,
                            int N_data_bins, G4double data_min_value, G4double data_max_value);
public:
    int _N_cols = 1;
    G4double col_min_value;
    G4double col_max_value;
    G4double col_value_width;

    int N_rows = 1;
    G4double row_min_value;
    G4double row_max_value;
    G4double row_value_width;

    int N_data_bins = 1;
    G4double data_min_value;
    G4double data_max_value;
    G4double data_bin_width;

    int N_total;
    G4double_vector histogram_matrix;
    //inline TwoDimensionalHistogram* getInstance(){return instance;};

//    void init_size(int N_cols, G4double col_min_value, G4double col_max_value,
//                    int N_rows, G4double row_min_value, G4double row_max_value,
//                    int N_data_bins, G4double data_min_value,G4double data_max_value);



    void accumulate(G4double row_value, G4double col_value, G4double data_value);

    inline void accumulate(int row_bin_Num, int col_bin_Num, int data_bin_Num);

    void accumulate(int direct_coordinate_bin, G4double data_value);

    inline void accumulate(int direct_coordinate_bin, int data_bin_Num);

//private:
    bool isInRange(G4double value, G4double min_value, G4double max_value);

    inline int getBinNumber(G4double value, G4double min_value, G4double bin_width)
    {
        return (int)((value-min_value)/bin_width);
    };


    inline int getIndex(int col_bin_num,int row_bin_num,int data_bin_num){
        return col_bin_num + (row_bin_num  + data_bin_num * N_rows) * _N_cols;
    };

    inline int getIndex4Data(int direct_coordinate_bin, int data_bin_num){
        return direct_coordinate_bin + data_bin_num * N_rows * _N_cols;
    };



};

#endif // TwoDimensionalHistogram_HH
