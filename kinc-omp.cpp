#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <vector>



struct PairwiseIndex
{
    int x;
    int y;
};



PairwiseIndex pairwise_index(size_t index)
{
    // compute pairwise index from scalar index
    size_t pos {0};
    size_t x {0};

    while ( pos + x <= index )
    {
        pos += x;
        ++x;
    }

    // return pairwise index
    return {
        static_cast<int>(x),
        static_cast<int>(index - pos)
    };
}



void operator++(PairwiseIndex& index)
{
    if ( ++index.y >= index.x )
    {
        index.y = 0;
        ++index.x;
    }
}



struct Pair
{
    PairwiseIndex index;
    float correlation;
};



void load_dataframe(const char *filename, int *p_rows, int *p_cols, float **p_data)
{
    // create input stream
    std::ifstream in(filename);

    // read sample names from first line
    std::string line;
    std::getline(in, line);

    // determine number of samples
    int cols {0};
    std::stringstream ss(line);

    // ignore RowID token
    std::string colname;
    ss >> colname;

    while ( !ss.eof() )
    {
        std::string colname;
        ss >> colname;

        cols++;
    }

    // read data from input file
    int rows {0};
    std::vector<float> values;

    while ( !in.eof() )
    {
        // read a line from the input file
        std::getline(in, line);

        std::stringstream ss(line);

        // read row name
        std::string rowname;
        ss >> rowname;

        // return if line is empty
        if ( ss.eof() )
        {
            break;
        }

        // read data elements
        for ( int i = 0; i < cols; i++ )
        {
            std::string token;
            ss >> token;

            // if token matches nan token then it as such
            if ( token == "NA" )
            {
                values.push_back(NAN);
            }

            // else this is a normal floating point expression
            else
            {
                float value;
                ss >> value;

                values.push_back(value);
            }
        }

        // increment number of rows
        rows++;
    }

    // initialize dataframe
    float *data = new float[rows * cols];

    memcpy(data, values.data(), rows * cols * sizeof(float));

    // save outputs
    *p_rows = rows;
    *p_cols = cols;
    *p_data = data;
}




float Similarity_compute(
    const float *        expressions,
    int                  n_samples,
    const PairwiseIndex& index,
    float                min_expression,
    float                max_expression,
    int                  min_samples)
{
    // initialize variables
    const float *x = &expressions[index.x * n_samples];
    const float *y = &expressions[index.y * n_samples];

    // compute intermediate sums
    int n = 0;
    float sumx = 0;
    float sumy = 0;
    float sumx2 = 0;
    float sumy2 = 0;
    float sumxy = 0;

    for ( int i = 0; i < n_samples; ++i )
    {
        // get pairwise sample
        float x_i = x[i];
        float y_i = y[i];

        // determine whether the sample is valid
        bool valid = 
            !std::isnan(x_i) &&
            !std::isnan(y_i) &&
            min_expression < x_i && x_i < max_expression &&
            min_expression < y_i && y_i < max_expression;

        // include pairwise sample if it is valid
        if ( valid )
        {
            sumx += x_i;
            sumy += y_i;
            sumx2 += x_i * x_i;
            sumy2 += y_i * y_i;
            sumxy += x_i * y_i;

            ++n;
        }
    }

    // compute correlation only if there are enough samples
    if ( n >= min_samples )
    {
        return (n*sumxy - sumx*sumy) / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
    }

    return NAN;
}



int main(int argc, char **argv)
{
    // parse command-line arguments
    if ( argc != 5 )
    {
        std::cerr << "usage: ./kinc-omp <bsize> <np> <infile> <outfile>\n";
        exit(-1);
    }

    int block_size = atoi(argv[1]);
    int np = atoi(argv[2]);
    const char *infile = argv[3];
    const char *outfile = argv[4];

    // set the number of openmp threads
    omp_set_num_threads(np);

    // measure execution time
    auto t1 = std::chrono::system_clock::now();

    // load dataframe
    int n_genes;
    int n_samples;
    float *expressions;

    load_dataframe(infile, &n_genes, &n_samples, &expressions);

    // printf("loaded dataframe (%d x %d)\n", n_genes, n_samples);

    // initialize execution parameters
    int min_samples {30};
    float min_expression {0};
    float max_expression {20};

    // initialize output
    std::vector<Pair> pairs;

    // iterate through all pairs
    int n_total_pairs = n_genes * (n_genes - 1) / 2;

    #pragma omp parallel for schedule(dynamic, 1)
    for ( int i = 0; i < n_total_pairs; i += block_size )
    {
        // initialize pairwise index
        PairwiseIndex index = pairwise_index(i);

        // process partition of pairs
        int n_pairs {std::min(block_size, n_total_pairs - i)};

        for ( int j = 0; j < n_pairs; ++j, ++index )
        {
            // execute similiarity kernel
            float correlation = Similarity_compute(
                expressions,
                n_samples,
                index,
                min_expression,
                max_expression,
                min_samples
            );

            // save results
            if ( !std::isnan(correlation) )
            {
                Pair pair { index, correlation };

                #pragma omp critical
                pairs.push_back(pair);
            }
        }
    }

    // save results
    std::ofstream out(outfile);

    for ( const Pair& p : pairs )
    {
        out
            << p.index.x << "\t"
            << p.index.y << "\t"
            << p.correlation << "\n";
    }
    
    // measure elapsed time
    auto t2 = std::chrono::system_clock::now();
    float elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0f;

    // print performance results
    std::cout
        << "kinc-omp" << "\t"
        << block_size << "\t"
        << np << "\t"
        << elapsed << "\n";

    return 0;
}
