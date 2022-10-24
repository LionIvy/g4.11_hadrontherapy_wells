#include "TwoDimensional_Histogram.hh"

TwoDimensional_Histogram::TwoDimensional_Histogram(int N_cols,      G4double col_min_value , G4double col_max_value,
                          int N_rows,      G4double row_min_value , G4double row_max_value,
                          int N_data_bins, G4double data_min_value, G4double data_max_value)
{
    G4double temp = 0;

    this->N_cols = N_cols;
    this->col_min_value = col_min_value;
    this->col_max_value = col_max_value;
    this->col_value_width = (col_max_value - col_min_value)/N_cols;
    temp = col_min_value;
    while(temp < col_max_value){
        columns_map.push_back(temp);
        temp += col_value_width;
    }

    this->N_rows = N_rows;
    this->row_min_value = row_min_value;
    this->row_max_value = row_max_value;
    this->row_value_width = (row_max_value - row_min_value)/N_rows;
    temp = row_min_value;
    while(temp < row_max_value){
        rows_map.push_back(temp);
        temp += row_value_width;
    }

    this->N_data_bins = N_data_bins;
    this->data_min_value = data_min_value;
    this->data_max_value = data_max_value;
    this->data_bin_width = (data_max_value - data_min_value)/N_data_bins;
    temp = data_min_value;
    while(temp < data_max_value){
        data_map.push_back(temp);
        temp += data_bin_width;
    }

    N_total = N_cols * N_rows * N_data_bins;

    G4double_vector null_vector(N_total, 0);
    histogram_matrix = null_vector;

};
//TwoDimensional_Histogram::~TwoDimensional_Histogram() {};


//void TwoDimensional_Histogram::init_size(int N_cols, G4double col_min_value, G4double col_max_value,
//                    int N_rows, G4double row_min_value, G4double row_max_value,
//                    int N_data_bins, G4double data_min_value,G4double data_max_value)
//{
//        this->N_cols = N_cols;
//        this->col_min_value = col_min_value;
//        this->col_max_value = col_max_value;
//        this->col_value_width = (col_max_value - col_min_value)/N_cols;

//        this->N_rows = N_rows;
//        this->row_min_value = row_min_value;
//        this->row_max_value = row_max_value;
//        this->row_value_width = (row_max_value - row_min_value)/N_rows;

//        this->N_data_bins = N_data_bins;
//        this->data_min_value = data_min_value;
//        this->data_max_value = data_max_value;
//        this->data_bin_width = (data_max_value - data_min_value)/N_data_bins;

//        N_total = N_cols * N_rows * N_data_bins;
//    };


void TwoDimensional_Histogram::accumulate(G4double row_value, G4double col_value, G4double data_value)
{
    if(isInRange(row_value,  row_min_value, row_max_value, row_value_width) &&
       isInRange(col_value,  col_min_value, col_max_value, col_value_width) &&
       isInRange(data_value,data_min_value,data_max_value,  data_bin_width))
    {
       histogram_matrix[getIndex(getBinNumber(row_value,  row_min_value, row_max_value, row_value_width),
                                 getBinNumber(col_value,  col_min_value, col_max_value, col_value_width),
                                 getBinNumber(data_value,data_min_value,data_max_value,  data_bin_width))] += 1;

    };
 };

void TwoDimensional_Histogram::accumulate(int row_bin_Num, G4double col_bin_Num, G4double data_bin_Num)
{

       histogram_matrix[getIndex(row_bin_Num,
                                 col_bin_Num,
                                 data_bin_Num)] += 1;
 };

bool TwoDimensional_Histogram::isInRange(G4double value, G4double min_value, G4double max_value, G4double value_step){
    if (value >= max_value){return false;}
    else if (value <  min_value){return false;}
    else{return true;};
}

int TwoDimensional_Histogram::getBinNumber(G4double value, G4double min_value, G4double max_value, G4double value_step){
    return (int)((value-min_value)/value_step);}
