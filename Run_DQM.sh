cmsRun test_cfg.py >> /dev/null
root -l -b -q DumpPlots_Digis.C
cp -r Digis /var/www/DQM/php-plots 
