
ADD_LIBRARY ( imathnumpymodule ${LIB_TYPE}
    imathnumpymodule.cpp
)

SET_TARGET_PROPERTIES ( imathnumpymodule
	PROPERTIES PREFIX "" SUFFIX ".so" BUILD_WITH_INSTALL_RPATH ON
)

INCLUDE_DIRECTORIES (
	${NUMPY_INCLUDE_DIRS}
)

TARGET_LINK_LIBRARIES ( imathnumpymodule
	PyImath
	PyIex
	Imath${ILMBASE_LIBSUFFIX}
	Iex${ILMBASE_LIBSUFFIX}
	${Boost_LIBRARIES}
)

INSTALL ( TARGETS imathnumpymodule
	DESTINATION lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages
)

