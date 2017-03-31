# this script is meant to be run from a machine on UW campus

output_dir=/astro/store/pogo4/danielsf/mlt_flares/

python assign_varParamStr.py --table stars_mlt_part_0870 \
--out_dir ${output_dir} --chunk_size 100000 --seed 88 &

sleep 5

python assign_varParamStr.py --table stars_mlt_part_1100 \
--out_dir ${output_dir} --chunk_size 100000 --seed 112 &

sleep 5

python assign_varParamStr.py --table stars_mlt_part_1160 \
--out_dir ${output_dir} --chunk_size 100000 --seed 321 &

sleep 5

python assign_varParamStr.py --table stars_mlt_part_1180 \
--out_dir ${output_dir} --chunk_size 100000 --seed 425 &

sleep 5

python assign_varParamStr.py --table stars_mlt_part_1220 \
--out_dir ${output_dir} --chunk_size 100000 --seed 6782 &

sleep 5

python assign_varParamStr.py --table stars_mlt_part_1250 \
--out_dir ${output_dir} --chunk_size 100000 --seed 371 &

sleep 5

python assign_varParamStr.py --table stars_mlt_part_1400 \
--out_dir ${output_dir} --chunk_size 100000 --seed 775 &

