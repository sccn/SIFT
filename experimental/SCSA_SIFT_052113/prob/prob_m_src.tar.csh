#!/bin/csh
#
#  Purpose:
#
#    Create a GZIP'ed TAR file of the m_src/prob files.
#
#  Modified:
#
#    02 January 2006
#
#  Author:
#
#    John Burkardt
#
#  Move to the directory just above the "prob" directory.
#
cd $HOME/public_html/m_src
#
#  Delete any TAR or GZ file in the directory.
#
echo "Remove TAR and GZ files."
rm prob/*.tar
rm prob/*.gz
#
#  Create a TAR file of the "prob" directory.
#
echo "Create TAR file."
tar cvf prob_m_src.tar prob/*
#
#  Compress the file.
#
echo "Compress the TAR file."
gzip prob_m_src.tar
#
#  Move the compressed file into the "prob" directory.
#
echo "Move the compressed file into the directory."
mv prob_m_src.tar.gz prob
#
#  Say goodnight.
#
echo "The prob_m_src gzip file has been created."
