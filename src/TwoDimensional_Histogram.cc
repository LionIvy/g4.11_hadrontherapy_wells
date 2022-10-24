#include "TwoDimensional_Histogram.hh"

//======================================================================================================================================
//                        TwoDimensionalHistogram
//======================================================================================================================================
TwoDimensionalHistogram::TwoDimensionalHistogram(int new_N_cols,      G4double new_col_min_value , G4double new_col_max_value,
                                                   int new_N_rows,      G4double new_row_min_value , G4double new_row_max_value,
                                                   int new_N_data_bins, G4double new_data_min_value, G4double new_data_max_value)
{
    init( new_N_cols,      new_col_min_value , new_col_max_value,
                        new_N_rows,      new_row_min_value , new_row_max_value,
                        new_N_data_bins, new_data_min_value, new_data_max_value);
}


//void TwoDimensionalHistogram::set_and_initialize(int new_N_cols,      G4double new_col_min_value , G4double new_col_max_value,
//                                                  int new_N_rows,      G4double new_row_min_value , G4double new_row_max_value,
//                                                  int new_N_data_bins, G4double new_data_min_value, G4double new_data_max_value)
//{
//    G4cout << "set_and_initialize cols N = " << this->N_cols << G4endl;
//    N_cols = new_N_cols;
//    G4cout << "set_and_initialize cols min" << G4endl;
//    col_min_value = new_col_min_value;
//    G4cout << "set_and_initialize cols max" << G4endl;
//    col_max_value = new_col_max_value;
//    G4cout << "set_and_initialize cols width" << G4endl;
//    col_value_width = (col_max_value - col_min_value)/N_cols;

//    G4cout << "set_and_initialize rows" << G4endl;
//    this->N_rows = new_N_rows;
//    this->row_min_value = new_row_min_value;
//    this->row_max_value = new_row_max_value;
//    this->row_value_width = (row_max_value - row_min_value)/N_rows;

//    G4cout << "set_and_initialize data" << G4endl;
//    this->N_data_bins = new_N_data_bins;
//    this->data_min_value = new_data_min_value;
//    this->data_max_value = new_data_max_value;
//    this->data_bin_width = (data_max_value - data_min_value)/N_data_bins;


//    N_total = N_cols * N_rows * N_data_bins;
//    G4cout << "null_vector(N_total, 0)" << G4endl;

//    //G4double_vector null_vector(N_total, 0);
//    G4cout << "histogram_matrix = null_vector;" << G4endl;
//    //histogram_matrix = null_vector;
//}

void TwoDimensionalHistogram::init(int N_cols_,      G4double _col_min_value , G4double _col_max_value,
                                   int _N_rows,      G4double _row_min_value , G4double _row_max_value,
                                   int _N_data_bins, G4double _data_min_value, G4double _data_max_value)
{
    G4cout << "set_and_initialize cols N = " << N_cols_ << G4endl;
    std::cout<<this->_N_cols<<std::endl;
    _N_cols = N_cols_;
    G4cout << "set_and_initialize cols min" << G4endl;
    col_min_value = _col_min_value;
    G4cout << "set_and_initialize cols max" << G4endl;
    col_max_value = _col_max_value;
    G4cout << "set_and_initialize cols width" << G4endl;
    col_value_width = (col_max_value - col_min_value)/_N_cols;

    G4cout << "set_and_initialize rows" << G4endl;
    this->N_rows = _N_rows;
    this->row_min_value = _row_min_value;
    this->row_max_value = _row_max_value;
    this->row_value_width = (row_max_value - row_min_value)/N_rows;

    G4cout << "set_and_initialize data" << G4endl;
    this->N_data_bins = _N_data_bins;
    this->data_min_value = _data_min_value;
    this->data_max_value = _data_max_value;
    this->data_bin_width = (data_max_value - data_min_value)/N_data_bins;


    N_total = _N_cols * N_rows * N_data_bins;
    G4cout << "null_vector(N_total, 0)" << G4endl;

    //G4double_vector null_vector(N_total, 0);
    G4cout << "histogram_matrix = null_vector;" << G4endl;
    //histogram_matrix = null_vector;
}

bool TwoDimensionalHistogram::isInRange(G4double value, G4double min_value, G4double max_value){
    if (value >= max_value){return false;}
    else if (value <  min_value){return false;}
    else{return true;};
}

void TwoDimensionalHistogram::accumulate(G4double row_value, G4double col_value, G4double data_value){
    if(isInRange(row_value,  row_min_value, row_max_value) &&
       isInRange(col_value,  col_min_value, col_max_value) &&
       isInRange(data_value,data_min_value,data_max_value)   )
    {
       histogram_matrix[getIndex(getBinNumber(row_value,  row_min_value, row_value_width),
                                 getBinNumber(col_value,  col_min_value, col_value_width),
                                 getBinNumber(data_value,data_min_value, data_bin_width))] += 1;

    };
}

inline void TwoDimensionalHistogram::accumulate(int row_bin_Num, int col_bin_Num, int data_bin_Num)
{

       histogram_matrix[getIndex(row_bin_Num,
                                 col_bin_Num,
                                 data_bin_Num)] += 1;
 }

void TwoDimensionalHistogram::accumulate(int direct_coordinate_bin, G4double data_value){
    if(isInRange(data_value,data_min_value,data_max_value))
    {
       histogram_matrix[getIndex4Data(direct_coordinate_bin,
                                      getBinNumber(data_value,data_min_value,  data_bin_width))] += 1;
    };
}

inline void TwoDimensionalHistogram::accumulate(int direct_coordinate_bin, int data_bin_num){
    histogram_matrix[getIndex4Data(direct_coordinate_bin, data_bin_num)] += 1;
}

