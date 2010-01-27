
CONFIGURE_FILE(doxygen.conf.in ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf)
ADD_CUSTOM_TARGET(docs
		      ${DOXYGEN} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf
          COMMAND cp ${CMAKE_SOURCE_DIR}/utils/doc/input/* ${CMAKE_CURRENT_BINARY_DIR}/doc/latex/
          COMMAND cp ${CMAKE_SOURCE_DIR}/utils/doc/input/* ${CMAKE_CURRENT_BINARY_DIR}/doc/
		      COMMAND make -C ${CMAKE_CURRENT_BINARY_DIR}/doc/latex
		      COMMAND mkdir -p ${CMAKE_SOURCE_DIR}/doc/html
		      COMMAND cp -r ${CMAKE_CURRENT_BINARY_DIR}/doc/html/* ${CMAKE_SOURCE_DIR}/doc/html/
 		      COMMAND cp    ${CMAKE_CURRENT_BINARY_DIR}/doc/latex/refman.pdf ${CMAKE_SOURCE_DIR}/doc/
 		      COMMAND cp    ${CMAKE_CURRENT_BINARY_DIR}/doc/latex/git_dales.pdf ${CMAKE_SOURCE_DIR}/doc/
 		      COMMAND cp    ${CMAKE_CURRENT_BINARY_DIR}/doc/latex/dales-manual.pdf ${CMAKE_SOURCE_DIR}/doc/
 		      COMMAND cp    ${CMAKE_CURRENT_BINARY_DIR}/doc/latex/dales-article.pdf ${CMAKE_SOURCE_DIR}/doc/
		      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf)
