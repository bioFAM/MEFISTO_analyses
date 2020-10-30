wget https://data.qiime2.org/distro/core/qiime2-2020.8-py36-linux-conda.yml
conda env create -n qiime2-2020.8 --file qiime2-2020.8-py36-linux-conda.yml
# OPTIONAL CLEANUP
rm qiime2-2020.8-py36-linux-conda.yml