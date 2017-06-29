# Rayden

[![Build status](https://ci.appveyor.com/api/projects/status/3sfvb2l0xe0x50ah?svg=true)](https://ci.appveyor.com/project/yakir12/rayden-jl) [![Build Status](https://travis-ci.org/yakir12/Rayden.jl.svg?branch=master)](https://travis-ci.org/yakir12/Rayden.jl)

[![Coverage Status](https://coveralls.io/repos/github/yakir12/Rayden.jl/badge.svg?branch=master)](https://coveralls.io/github/yakir12/Rayden.jl?branch=master) [![codecov.io](http://codecov.io/github/yakir12/Rayden.jl/coverage.svg?branch=master)](http://codecov.io/github/yakir12/Rayden.jl?branch=master)

This Julia package allows for geometric ray tracing with ellipsoids. It includes intersection and refraction/reflection of rays with arbitrary ellipsoids. It accomplishes that in about 100 lines of code thanks to heavy use of `CoordinateTransformations.jl` and `StaticArrays.jl`.
