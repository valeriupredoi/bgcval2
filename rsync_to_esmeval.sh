# Script to copy report to esmeval public facing server.
# You may need to request access to the esmeval CEDA-JASMIN group
# or set your own path.
# old: /gws/nopw/j04/esmeval/public/CompareReports/bgcval2/
Outpath=/gws/ssde/j25a/esmeval/public/CompareReports/bgcval2/$USER
mkdir -p $Outpath
rsync -av CompareReports2/* $Outpath/.

# this outputs to the web address:
# https://gws-access.jasmin.ac.uk/public/esmeval/CompareReports/bgcval2/$USER
