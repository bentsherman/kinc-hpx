#include <hpx/hpx_init.hpp>
#include <hpx/hpx.hpp>

#include <hpx/include/parallel_algorithm.hpp>
#include <boost/range/irange.hpp>

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <utility>
#include <sstream>
#include <string>
#include <vector>



struct PairwiseIndex
{
public:
    PairwiseIndex() {}

    PairwiseIndex(size_t index)
    {
        // compute pairwise index from scalar index
        size_t pos {0};
        size_t x {0};

        while ( pos + x <= index )
        {
            pos += x;
            ++x;
        }

        this->x = static_cast<int>(x);
        this->y = static_cast<int>(index - pos);
    }

public:
    int x;
    int y;
};



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



void load_dataframe(const std::string& filename, int *p_rows, int *p_cols, float **p_data)
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



struct Block
{
public:
    Block(int start, int size):
        _start(start),
        _size(size)
    {}

    int start() const { return _start; }
    int size() const { return _size; }
    const std::vector<Pair>& pairs() const { return _pairs; }
    std::vector<Pair>& pairs() { return _pairs; }

private:
    int _start;
    int _size;
    std::vector<Pair> _pairs;
};



struct TaskRunner
{
    typedef hpx::future<Block> BlockFuture;

    static Block process(
        const Block & block,
        const float * expressions,
        int           n_samples,
        float         min_expression,
        float         max_expression,
        int           min_samples)
    {
        // initialize result block
        Block result(block.start(), block.size());

        // initialize pairwise index
        PairwiseIndex index(block.start());

        // process partition of pairs
        for ( int i = 0; i < block.size(); ++i, ++index )
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
                result.pairs().push_back({
                    index,
                    correlation
                });
            }
        }
        
        return result;
    }

    std::vector<Pair> execute(
        const float * expressions,
        int           n_genes,
        int           n_samples,
        float         min_expression,
        float         max_expression,
        int           min_samples,
        int           block_size)
    {
        std::vector<BlockFuture> futures;

        // initiate work blocks
        int n_total_pairs = n_genes * (n_genes - 1) / 2;

        for ( int i = 0; i < n_total_pairs; i += block_size )
        {
            futures.push_back(hpx::async(
                process,
                Block(i, block_size),
                expressions,
                n_samples,
                min_expression,
                max_expression,
                min_samples
            ));
        }

        // aggregate results
        std::vector<Pair> pairs;

        for ( auto& block_future : futures )
        {
            auto block {block_future.get()};
            pairs.insert(pairs.end(), block.pairs().begin(), block.pairs().end());
        }

        return pairs;
    }
};



int hpx_main(hpx::program_options::variables_map& vm)
{
    // initialize execution parameters
    std::string infile = vm["infile"].as<std::string>();
    std::string outfile = vm["outfile"].as<std::string>();
    int block_size = vm["bsize"].as<int>();
    int min_samples {30};
    float min_expression {0};
    float max_expression {20};

    // initialize task runner
    TaskRunner runner;

    // measure execution time
    std::uint64_t t = hpx::util::high_resolution_clock::now();

    // load dataframe
    int n_genes;
    int n_samples;
    float *expressions;

    load_dataframe(infile, &n_genes, &n_samples, &expressions);

    // printf("loaded dataframe (%d x %d)\n", n_genes, n_samples);

    // execute similarity
    auto pairs = runner.execute(
        expressions,
        n_genes,
        n_samples,
        min_expression,
        max_expression,
        min_samples,
        block_size
    );

    // save results
    std::ofstream out(outfile);

    for ( const Pair& p : pairs )
    {
        out
            << p.index.x << "\t"
            << p.index.y << "\t"
            << p.correlation << "\n";
    }

    double elapsed = (hpx::util::high_resolution_clock::now() - t) / 1e9;

    // print performance results
    std::uint64_t np = hpx::get_os_thread_count();

    std::cout
        << "kinc-hpx" << "\t"
        << block_size << "\t"
        << np << "\t"
        << elapsed << "\n";

    return hpx::finalize();
}



int main(int argc, char* argv[])
{
    using namespace hpx::program_options;

    // Configure application-specific options.
    options_description desc_commandline;

    desc_commandline.add_options()
        ("infile", value<std::string>(), "Input expression matrix")
        ("outfile", value<std::string>(), "Output correlation matrix")
        ("bsize", value<int>()->default_value(1024), "Work block size")
    ;

    // Initialize and run HPX
    return hpx::init(desc_commandline, argc, argv);
}
