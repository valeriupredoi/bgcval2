import iris
import numpy as np
import iris.coord_systems
import iris.fileformats
from cf_units import Unit


fx_data = np.empty((4, 3, 3))
fx_data_3 = np.empty((4, 3, 3, 3))
fx_data_2 = np.empty((3, 3))
fx_data[:] = 60.
fx_data[1, 2] = 30.
new_cube_data = np.empty((2, 3, 3))
new_cube_data[:] = 200.
crd_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
lons = iris.coords.DimCoord([0, 10.5, 219],
                            # standard_name='longitude',
                            var_name='lon',
                            bounds=[[0, 1], [10, 20], [210, 220]],
                            units='degrees_east',
                            coord_system=crd_sys)
lats = iris.coords.DimCoord([0, 70.5, 130],
                            # standard_name='latitude',
                            var_name='lat',
                            bounds=[[0, 1], [70, 71], [120, 130]],
                            units='degrees_north',
                            coord_system=crd_sys)
levs = iris.coords.DimCoord([0, 10.5, 30],
                            var_name='depth',
                            bounds=[[0, 1], [10, 20], [20, 30]],
                            units='m')
times = iris.coords.DimCoord([0, 1.5, 2.5, 3.5],
                             standard_name='time',
                             bounds=[[0, 1], [1, 2], [2, 3],
                                     [3, 4]],
                             units=Unit('days since 1950-01-01 00:00:00', calendar='standard'))
fx_coords_spec = [(times, 0), (lats, 1), (lons, 2)]
fx_coords_spec_3 = [(times, 0), (levs, 1), (lats, 2), (lons, 3)]
fx_coords_spec_2 = [(lats, 0), (lons, 1)]
fx_mask = iris.cube.Cube(fx_data,
                         dim_coords_and_dims=fx_coords_spec,
                         var_name="tmask",
                         units='%')
fx_mask_2 = iris.cube.Cube(fx_data_2,
                         dim_coords_and_dims=fx_coords_spec_2,
                         var_name="area",
                         units='%')
e1t = iris.cube.Cube(fx_data,
                         dim_coords_and_dims=fx_coords_spec,
                         var_name="e1t",
                         units='%')
e2t = iris.cube.Cube(fx_data,
                         dim_coords_and_dims=fx_coords_spec,
                         var_name="e2t",
                         units='%')
e2u = iris.cube.Cube(fx_data,
                         dim_coords_and_dims=fx_coords_spec,
                         var_name="e2u",
                         units='%')
e3v = iris.cube.Cube(fx_data,
                         dim_coords_and_dims=fx_coords_spec,
                         var_name="e3v",
                         units='%')
e3u = iris.cube.Cube(fx_data_3,
                         dim_coords_and_dims=fx_coords_spec_3,
                         var_name="e3u",
                         units='%')
u3d = iris.cube.Cube(fx_data_3,
                         dim_coords_and_dims=fx_coords_spec_3,
                         var_name="u3d",
                         units='%')
e1v = iris.cube.Cube(fx_data,
                         dim_coords_and_dims=fx_coords_spec,
                         var_name="e1v",
                         units='%')
umask = iris.cube.Cube(fx_data,
                         dim_coords_and_dims=fx_coords_spec,
                         var_name="umask",
                         units='%')
navlat = iris.cube.Cube([0, 1.5, 3],
                         dim_coords_and_dims=(),
                         var_name="nav_lat",
                         units='degrees')
soice = iris.cube.Cube(fx_data,
                         dim_coords_and_dims=fx_coords_spec,
                         var_name="soicecov",
                         units='%')
iris.save(iris.cube.CubeList([fx_mask, fx_mask_2, e1t, e2t, navlat, e2u, e3v, e1v, umask]), "local_test/mesh_mask_eORCA1_wrk.nc")
iris.save(iris.cube.CubeList([fx_mask,soice,e3u,u3d]), "local_test/BGC_data/u-cp416debug/nemo_u-cp416debugo_1y_1990_grid-T.nc")
iris.save(iris.cube.CubeList([fx_mask,soice,e3u,u3d]), "local_test/BGC_data/u-cp416debug/nemo_u-cp416debugo_1y_1990_grid-U.nc")
iris.save(iris.cube.CubeList([fx_mask,soice,e3u,u3d]), "local_test/BGC_data/u-cp416debug/nemo_u-cp416debugo_1y_1990_grid-V.nc")
iris.save(fx_mask, "local_test/bgc/WOA/annual/woa13_decav_t00_01v2.nc")
iris.save(fx_mask, "local_test/bgc/WOA/annual/woa13_decav_s00_01v2.nc")
