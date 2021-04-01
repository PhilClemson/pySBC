# pySBC
A Python implementation of simulation based calibration for CmdStan.

# Usage
Place pySBC.py in the main directory of CmdStan and run on the command line using `python pySBC.py --option []`.

## Options
`--exe` Location of compiled Stan executable for computing SBC. The rank statistics should be defined in an array named `lt_sim` in the `generated quantities` block.

`--nt` Number of threads.

`--warmup` Number of warmup samples.

`--samples` Number of samples after warmup.

`--data` Location of data file (use `""` for no data file).

`--J` Number of bins in the rank statistic histogram.

`--ej` Expected number in each rank statistic bin (if uniformly distributed). Should be no less than 5.

`--output` Name of output file to write histogram data to.

`--t` Maximum time to wait for a chain to complete in seconds.

## Example
`python pySBC.py --exe /benchmarks/eight_schools/eight_schools_sbc --nt 8 --warmup 1000 --samples 1000 --data /benchmarks/eight_schools/eight_schools.data.R --J 10 --ej 5 --output sbc_eight_schools_histogram.out --t 10`
