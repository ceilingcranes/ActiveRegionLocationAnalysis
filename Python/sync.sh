#!/bin/bash
# Small script to sync Python code to GPU machine.
user=maha7656
rsync -rP ./*.py $user@swami.colorado.edu:/home/$user/thesis_work/code/Python --exclude=.env
# rsync -rP data/distance_data.csv $user@swami.colorado.edu:/home/$user/thesis_work/code/Python --exclude=.env
