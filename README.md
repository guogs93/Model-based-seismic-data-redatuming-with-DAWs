# Redatuming with DAWs

This repository contains the source code for model-based seismic data redatuming using time-domain data-assimilated wavefield reconstruction (DAWs).

The code was developed by **Gaoshan Guo** during his PhD (2021.01~2024.03) at the **CNRS–Géoazur laboratory**, under the supervision of **Stéphane Operto**, and supported by the **WIND consortium** ([https://www.geoazur.fr/WIND/bin/view](https://www.geoazur.fr/WIND/bin/view)).

We gratefully acknowledge the sponsorship of **AkerBP, ExxonMobil, Petrobras, Shell**, and **SINOPEC**.

The method implements advanced time-domain seismic redatuming with data-assimilated wavefield reconstruction. It is designed for high-performance computing environments and has been tested on various supercomputing platforms.

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
