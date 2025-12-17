# MetaStrain

**Efficient Cell Factory Design by Combining Meta-Heuristic Algorithm with Enzyme Constrained Metabolic Models**

## Overview

MetaStrain is a computational framework for optimal strain design that integrates meta-heuristic optimization algorithms with enzyme-constrained genome-scale metabolic models (ecGEMs). The framework employs the JADE (Adaptive Differential Evolution with Optional External Archive) algorithm to identify optimal gene manipulation strategies (overexpression, knockdown, and knockout) for maximizing target product yields in microbial cell factories.

## Key Features

- **Meta-heuristic Optimization**: Implements JADE algorithm for efficient exploration of the combinatorial gene manipulation space
- **Enzyme-Constrained Models**: Supports enzyme-constrained genome-scale metabolic models (ecGEMs) for more realistic predictions
- **Multiple Gene Operations**: Simultaneously considers overexpression, knockdown, and knockout strategies
- **Parallel Computing**: Leverages Ray framework for distributed fitness evaluation
- **Dimensionality Reduction**: Utilizes pre-processed target gene sets from ecFSEOF analysis to reduce search space
- **Multiple Organisms**: Supports both *E. coli* (ec_iML1515) and *S. cerevisiae* (ecYeast) models
- **MOMA Integration**: Incorporates Minimization of Metabolic Adjustment (MOMA) for flux prediction under genetic perturbations

## Project Structure

```
MetaStrain/
├── Cases/                          # Entry point scripts for different case studies
│   ├── start_ec_iML1515.py        # E. coli tryptophan production case
│   ├── start_JADE_artemisinic_acid.py  # Yeast artemisinic acid production case
│   ├── start_JADE_spermidine.py   # Yeast spermidine production case
│   └── ModelOperation_eciML1515.py    # Model operations for E. coli
│
├── models/                         # Metabolic models and pre-processed data
│   ├── dimensionality_reduction/   # Reduced target gene sets from ecFSEOF
│   │   ├── ecFSEOF_MAPPING_eciML1515*.csv
│   │   ├── ecFSEOF_MAPPING_AA.csv
│   │   └── ecFSEOF_MAPPING_sp.csv
│   ├── eciML1515_batch.mat        # E. coli enzyme-constrained model
│   ├── fixed_eciML1515_batch.mat
│   ├── ecYeastGEM_batch.mat       # Yeast enzyme-constrained model
│   └── *.xls, *.xml, *.json       # Additional model formats
│
├── search module/                  # Core optimization algorithms
│   ├── JADE_operator.py           # JADE algorithm implementation
│   ├── GA_Operator.py             # Genetic algorithm operators
│   ├── ModelOperation.py          # Model manipulation functions
│   ├── MPMA.py                    # MOMA implementation
│   └── start_JADE.py              # Generic JADE execution script
│
├── tools/                          # Analysis and utility tools
│   ├── FVAanalysis.py             # Flux variability analysis
│   ├── Get_flux.py                # Flux extraction utilities
│   ├── gene_filter.py             # Gene filtering functions
│   ├── mapping.py                 # Gene-protein mapping
│   ├── phaseplane.py              # Phase plane analysis
│   └── Single_target.py           # Single target analysis
│
├── Result_analysis/                # Result analysis scripts
│   ├── result_checking.py
│   └── result_count.py
│
└── Redundancy_analysis/            # Redundancy and essentiality analysis
    ├── essential_genes.py
    └── enforce_Ntarget_selector.py
```

## Installation

### Prerequisites

- Python 3.7 or higher
- Gurobi Optimizer (for linear programming solver)
- Required Python packages (see below)

### System Requirements

- **Memory**: Recommended 8GB RAM or more (depends on model size)
- **CPU**: Multi-core processor recommended for parallel execution
- **Storage**: ~500MB for models and dependencies
- **Operating System**: Linux, macOS, or Windows (tested on Linux)

### Dependencies

Install the required Python packages:

```bash
pip install cobra
pip install ray
pip install numpy
pip install pandas
pip install matplotlib
pip install gurobipy
pip install openpyxl  # For Excel file reading
```

**Note**: Gurobi requires a license. Academic licenses are available free of charge. Please refer to [Gurobi's website](https://www.gurobi.com/) for license installation instructions.

## Quick Start

### Running a Case Study

To run a specific case study, navigate to the `Cases` directory and execute the corresponding `start_*.py` script:

#### Example 1: E. coli Tryptophan Production

```bash
cd Cases
python start_ec_iML1515.py
```

This script will:
1. Load the ec_iML1515 model
2. Set up the optimization problem for tryptophan production
3. Run JADE algorithm to identify optimal gene manipulation strategies
4. Save results to CSV files (`JADE_ind_trp.csv` and `JADE_val_trp.csv`)
5. Display convergence curves

#### Example 2: Yeast Artemisinic Acid Production

```bash
cd Cases
python start_JADE_artemisinic_acid.py
```

#### Example 3: Yeast Spermidine Production

```bash
cd Cases
python start_JADE_spermidine.py
```

### Configuration

Key parameters can be adjusted in each `start_*.py` script:

- **`n`**: Number of target genes (dimension of optimization problem)
- **`popsize`**: Population size for JADE algorithm (default: 100)
- **`totalTime`**: Number of independent runs (default: 1-2)
- **`FES`**: Maximum function evaluations (default: 200-500)
- **`c`**: Learning rate for parameter adaptation (default: 0.1)
- **`p`**: Percentage of best individuals for p-best selection (default: 0.05)

## Methodology

### Workflow

1. **Model Preparation**: Load enzyme-constrained metabolic model and set growth/product constraints
2. **Target Gene Selection**: Load pre-processed target gene sets from `dimensionality_reduction/` folder (obtained via ecFSEOF analysis)
3. **Reference Solution**: Calculate wild-type flux distribution using MOMA
4. **Optimization**: Apply JADE algorithm to search for optimal gene manipulation strategies
5. **Evaluation**: For each candidate solution, modify enzyme constraints and evaluate product yield using MOMA
6. **Output**: Save best solutions and convergence history

### Gene Manipulation Encoding

Each gene can be assigned one of four discrete states:
- **0**: No adjustment (wild-type)
- **1**: Overexpression (increase enzyme lower bound to 4× reference flux)
- **2**: Knockdown (decrease enzyme upper bound to 0.5× reference flux)
- **3**: Knockout (set enzyme upper bound to 0)

The JADE algorithm operates on continuous values [0, 4] which are mapped to discrete operations during fitness evaluation.

### Fitness Function

The fitness is defined as the product yield (product flux / substrate uptake flux) calculated using MOMA:

```
fitness = product_flux / substrate_uptake_flux
```

The optimization maximizes this fitness value (minimizes the negative fitness).

## Output Files

Each run generates two main output files:

1. **`JADE_ind_*.csv`**: Contains the best individuals (gene manipulation strategies) found in each independent run, with columns:
   - `Gene_1`, `Gene_2`, ..., `Gene_n`: Gene manipulation states (0-3)
   - `Best Fitness`: Maximum fitness value achieved

2. **`JADE_val_*.csv`**: Contains the convergence history, with columns:
   - `Iteration_0`, `Iteration_1`, ...: Fitness values at each generation

## Model Files

### Supported Models

- **ec_iML1515**: Enzyme-constrained *E. coli* iML1515 model
- **ecYeastGEM**: Enzyme-constrained *S. cerevisiae* model

### Target Gene Sets

Pre-processed target gene sets are stored in `models/dimensionality_reduction/`:
- `ecFSEOF_MAPPING_eciML1515*.csv`: Target genes for E. coli tryptophan production
- `ecFSEOF_MAPPING_AA.csv`: Target genes for artemisinic acid production
- `ecFSEOF_MAPPING_sp.csv`: Target genes for spermidine production

These files contain gene names and corresponding enzyme reaction IDs, pre-filtered using ecFSEOF (enzyme-constrained Flux Scanning based on Enforced Objective Flux) analysis.

## Advanced Usage

### Customizing for New Products

To adapt MetaStrain for a new product:

1. **Prepare Model**: Ensure your enzyme-constrained model is in the correct format (MATLAB `.mat` file)
2. **Identify Target Genes**: Run ecFSEOF analysis to identify candidate target genes
3. **Create Case Script**: Copy an existing `start_*.py` script and modify:
   - Model loading function
   - Product reaction ID
   - Growth/product constraints
   - Target gene file path
   - Gene-protein mapping

### Parallel Execution

The framework uses Ray for parallel fitness evaluation. To adjust the number of workers:

```python
ray.init(num_cpus=N)  # Use N CPU cores
```

### Understanding Results

The optimization results represent gene manipulation strategies where:
- Each row in `JADE_ind_*.csv` corresponds to one independent run
- Gene values (0-3) indicate the manipulation type for each target gene
- Higher fitness values indicate better product yields
- Multiple runs help assess solution robustness and convergence

## Troubleshooting

### Common Issues

1. **Gurobi License Error**
   - Ensure Gurobi license is properly installed and activated
   - Check license file location: `~/.gurobi/gurobi.lic`
   - For academic use, obtain free license from Gurobi website

2. **Model Loading Errors**
   - Verify model file paths are correct (use relative paths from script location)
   - Ensure model files are in the `models/` directory
   - Check that model format matches expected type (`.mat` for MATLAB models)

3. **Ray Initialization Issues**
   - If Ray fails to initialize, try: `ray.init(ignore_reinit_error=True)`
   - Check available system resources (CPU cores, memory)
   - On some systems, may need to specify: `ray.init(address='local')`

4. **Infeasible Solutions**
   - Many infeasible solutions (fitness = -0.0001) may indicate:
     - Overly restrictive constraints
     - Incompatible gene manipulations
     - Model inconsistencies
   - Try relaxing growth/product constraints or adjusting gene manipulation bounds

5. **Path Issues**
   - Ensure script is run from correct directory
   - Use relative paths or update absolute paths in scripts
   - Check that all required files exist in expected locations

### Performance Optimization

- **Population Size**: Larger populations (100-200) improve exploration but increase computation time
- **Function Evaluations**: More evaluations (500-1000) may improve solution quality for complex problems
- **Parallel Workers**: Set `ray.init(num_cpus=N)` to match available CPU cores for optimal speedup
- **Model Tolerance**: Default tolerance (1e-9) balances accuracy and speed; can be adjusted if needed
- **Multiple Runs**: Running multiple independent runs (`totalTime > 1`) helps identify robust solutions

## Citation

If you use MetaStrain in your research, please cite:

```
[Authors]. MetaStrain: Efficient Cell Factory Design by Combining Meta-Heuristic Algorithm 
with Enzyme Constrained Metabolic Models. bioRxiv (2025). 
doi: https://doi.org/10.1101/2025.06.12.659423v2
```

Or in BibTeX format:

```bibtex
@article{metastrain2025,
  title={MetaStrain: Efficient Cell Factory Design by Combining Meta-Heuristic Algorithm with Enzyme Constrained Metabolic Models},
  author={[Authors]},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2025.06.12.659423v2}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**MIT License**

Copyright (c) 2025 [Your Name/Institution]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Contact

For questions, issues, or contributions, please contact: mbliao0824@163.com

## Acknowledgments

This work integrates:
- COBRApy for metabolic modeling
- Ray for distributed computing
- Gurobi for optimization
- GECKO for enzyme-constrained modeling

## References

1. Sánchez, B. J., et al. (2017). Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic constraints. *Molecular Systems Biology*, 13(8), 935.

2. Zhang, J., & Sanderson, A. C. (2009). JADE: Adaptive differential evolution with optional external archive. *IEEE Transactions on Evolutionary Computation*, 13(5), 945-958.

3. Segrè, D., et al. (2002). Analysis of optimality in natural and perturbed metabolic networks. *Proceedings of the National Academy of Sciences*, 99(23), 15112-15117.
......
