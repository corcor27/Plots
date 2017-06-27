# Input inference folder
FOLDER=$1

ITER=11999 #Take last iteration - change this if you run with more
INPUT_FILE=${FOLDER}/output.hdf
OUTPUT_FILE=${FOLDER}/posteriors.png
pycbc_inference_plot_posterior \
    --iteration ${ITER} \
    --input-file ${INPUT_FILE} \
    --output-file ${OUTPUT_FILE} \
    --plot-scatter \
    --plot-marginal \
    --z-arg logplr \
    --parameters mass1 mass2 \
                 mchirp q \
                 chi_p chi_eff \
