RUN=$1

python whitened.py --input-file final/${RUN}/output.hdf \
                   --plot-map-waveforms \
                   --min-xlim -2.5 \
                   --injection-file final/${RUN}/injection.xml.gz \
                   --output-file final/${RUN}/strain_plot25

python whitened.py --input-file final/${RUN}/output.hdf \
                   --plot-map-waveforms \
                   --min-xlim -1.5 \
                   --injection-file final/${RUN}/injection.xml.gz \
                   --output-file final/${RUN}/strain_plot15

python whitened.py --input-file final/${RUN}/output.hdf \
                   --plot-map-waveforms \
                   --min-xlim -1.0 \
                   --injection-file final/${RUN}/injection.xml.gz \
                   --output-file final/${RUN}/strain_plot10

python whitened.py --input-file final/${RUN}/output.hdf \
                   --plot-map-waveforms \
                   --min-xlim -0.5 \
                   --injection-file final/${RUN}/injection.xml.gz \
                   --output-file final/${RUN}/strain_plot05

python whitened.py --input-file final/${RUN}/output.hdf \
                   --plot-map-waveforms \
                   --min-xlim -0.2 \
                   --injection-file final/${RUN}/injection.xml.gz \
                   --output-file final/${RUN}/strain_plot02

# Rename dictionary
mv final/${RUN}/strain_plot_dic.npy final/${RUN}/mapDic.npy
