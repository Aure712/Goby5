bash liftOver.sh > 1_liftOver.log 2>&1
echo "liftOver.sh completed at $(date)"
bash convert.sh > 2_convert.log 2>&1
echo "liftOver.sh completed at $(date)"
