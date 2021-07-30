docker run -it --rm --gpus all -u $(id -u):$(id -g) \
    -v $(pwd):/run  \
    -v /home/huhlim/junk/alphafold:/home/huhlim/junk/alphafold \
    -v /feig/s1/huhlim/apps/AlphaFold:/feig/s1/huhlim/apps/AlphaFold \
    -v /home/huhlim/apps:/home/huhlim/apps \
    -v /home/huhlim/db:/home/huhlim/db \
    -v /feig/s1/huhlim/db:/feig/s1/huhlim/db \
    -e CUDA_VISIBLE_DEVICES \
    alphafold_debug
