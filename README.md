# Redatuming with DAWs

This repository contains the source code for model-based seismic data redatuming using time-domain data-assimilated wavefield reconstruction (DAWs).

This code was developed by Gaoshan Guo during the PHD study at CNRS-GEOAZUR labotary, under the supervision of Stephane Operto, supported by the WIND consortium (https://www.geoazur.fr/WIND/bin/view).

We thanks for the sponser of AkerBP, ExxonMobil, Petrobras, Shell, and SINOPEC.

The method implements advanced time-domain seismic redatuming with data-assimilated wavefield reconstruction. It is designed for high-performance computing environments and tested on various supercomputing platforms.

## Usage

This code is intended for researchers working in seismic redatuming of marine seismic data. Compilation and running instructions depend on the computing environment—please refer to your corresponding `Makefile.inc_*` file.

To clean compiled objects and binaries, run:
make clean && make all

## Citation

If you use this code in your research or publication, please cite the following paper:

Gaoshan Guo, Stéphane Operto, Model-based seismic data redatuming with time-domain data-assimilated wavefield reconstruction, Geophysical Journal International, 2025, in press.

## License
Academic Use Only
This software is provided for non-commercial academic research purposes only.
Commercial use is strictly prohibited without explicit permission from the authors.

See the LICENSE file for details.
