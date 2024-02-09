YEAR=2022a
VERSION=v1


XSECFILE=xsec_covariance_DUNE_systs_${YEAR}_${VERSION}_XsecOnly.xml
FLUXFILE=DUNE_flux_systs.xml
COMBFILE=xsec_covariance_DUNE_systs_${YEAR}_${VERSION}.xml

cp ${FLUXFILE} ${FLUXFILE}.tmp
sed -i '1,2d' ${FLUXFILE}.tmp

cp ${XSECFILE} ${COMBFILE}
sed -i '$d' ${COMBFILE}

echo $'\n' >> ${COMBFILE}
cat ${FLUXFILE}.tmp >> ${COMBFILE}
rm ${FLUXFILE}.tmp

echo "Made file:" ${COMBFILE}
