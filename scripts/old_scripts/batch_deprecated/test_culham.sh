# Deprecated script.

# Original Output:
# terminate called after throwing an instance of 'std::runtime_error'
# what():  Wrong choice for the geometry

if [ -z $1 ] || [ -z $2 ]
then
  echo "Please provide cell sizes for test"
else
  if [ $(($1 % 2)) -eq 1 ]
  then
    echo "There must be an even number of cells (r direction)"
  else
    if [ $(($2 % 2)) -eq 1 ]
    then
      echo "There must be an even number of cells (theta direction)"
    else

      mkdir -p culham

      cd culham

      python3 ../../python/make_points_r.py $1 r.txt --R0 1e-3 --R 1.0
      python3 ../../python/make_points_theta.py $2 theta.txt

      build_dir=build_gnu

      # Refactoring changes:
      # ../${build_dir}/main -n 4 -a 0 --mod_pk 3 --DirBC_Interior 1 --divideBy2 0 -r 1e-5  --smoother 3 --verbose 2 --debug 0 --extrapolation 0 --optimized 1 --openmp 4 --v1 1 --v2 1 -R 1 --prob 7  --maxiter 300 --alpha_coeff 2 --beta_coeff 1 --res_norm 3 --f_grid_r "r.txt" --f_grid_theta "theta.txt" --rel_red_conv 1e-9 --f_sol_out out.txt
      ../${build_dir}/gmgpolar --nr_exp 4 --anisotropic_factor 0 --geometry 3 --DirBC_Interior 1 --divideBy2 0 --R0 1e-5 --verbose 2 --extrapolation 0 --openmp 4 --preSmoothingSteps 1 --postSmoothingSteps 1 --Rmax 1.0 --problem 1 --maxIterations 300 --alpha_coeff 3 --alpha_jump 0.7081 --beta_coeff 1 --residualNormType 2 --write_grid_file $write_grid_file --load_grid_file 1 --file_grid_radii "r.txt" --file_grid_angles "theta.txt" 

      python3 ../../python/plot_culham.py

    fi
  fi
fi
