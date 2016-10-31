# Installation Guide #
1. Download pywikibot developper version 3.0-dev from  
http://tools.wmflabs.org/pywikibot/  
or use this direct link to download: tools.wmflabs.org/pywikibot/core_old.tar.gz
Note: annlightment won't work with versions 2.0 or lower
2. Extract core_old.tar.gz to any location. Go to this location (core_old)
and make a link of its subfolder pywikibot
3. Copy the link into annlightment folder and name the link 'pywikibot'
4. Make sure the annlightment folder contains a 'user-config.py' (It should already exist)
5. Now copy the family file 'TillsWiki_family.py' (located in annlightment folder) to core_old/pywikibot/families
6. To start the programm got the annlightment folder and run annlightment.py:
$ python3 annlightment.py upload path_to_ANNOgesic_operon_gff \  
path_to_ANNOgesic_merge.csv Q###  
You can also run python3 annlightment.py -h to get help. Annlightment may ask you to enter the passphrase for the Bot account.
