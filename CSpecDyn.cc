#include "CSpecDyn.h"

// Fyi: Immer FORTRAN order!!! (stride1 ist nicht definiert!)

CSpecDyn::CSpecDyn():
N(NUM), pdims(PDIMS), cfl(CFL), out_dir(OUT_DIR), out_interval(OUT_INTERVAL), end_simu(END_SIMU), L(LENGTH), nu(NU), eta(ETA), setup(SETUP)
{
  
  // init MPI
  MPI_Init(NULL, NULL);
 
  // init FFT  
  FFT = MikeFFT(N, pdims);
  
  // get sizes
  FFT.get_sizes(size_R, start_R, size_F, start_F);
  
  size_R_tot = size_R[0]*size_R[1]*size_R[2];
  size_F_tot = size_F[0]*size_F[1]*size_F[2];
  
  // get MPI info
  FFT.get_comm(comm);
  
  MPI_Comm_size(comm, &nprocs);
	MPI_Comm_rank(comm, &myRank);
  MPI_Cart_coords(comm, myRank, 2, mpi_coords);
  
  // drived quantities
  XB = -L*0.5;
  dx = L/N;
  dk = PI2/L;
  dt = cfl*dx;
  time = 0.;
  
  // allocate memory
  kx = new double[size_F[0]];
  ky = new double[size_F[1]];
  kz = new double[size_F[2]];
  k2 = new double[size_F_tot];
  
  Vx_R = FFT.malloc_R();
  Vy_R = FFT.malloc_R();
  Vz_R = FFT.malloc_R();
  
  Bx_R = FFT.malloc_R();
  By_R = FFT.malloc_R();
  Bz_R = FFT.malloc_R();
  
  float_array = (float*) malloc(sizeof(float)*size_R_tot);
  float_array_vector = (float*) malloc(sizeof(float)*size_R_tot*3);
  
  // setup initial fields
  setup_k();
  setup_fields();
  
  // create output directory
  if(myRank == 0)
	{
    mkdir(out_dir.c_str(), 0777);
	}MPI_Barrier(comm);
  
  // create subarray for mpi-vti output
  int size_total[3] = {N,N,N};
  MPI_Type_create_subarray(3, size_total, size_R, start_R, MPI_ORDER_C, MPI_FLOAT, &vti_subarray);
  MPI_Type_commit(&vti_subarray);
  
}

void CSpecDyn::setup_k()
{
  // kx
  for(int ix = 0; ix < size_F[0]; ix++)
	{
    
    int ix_glob = start_F[0] + ix;
    
    if(ix_glob < N/2)
    {
      kx[ix] = ix_glob * dk;
    }
    else
    {
      kx[ix] = (ix_glob-N) * dk;
		}
	}
  
  // ky
  for(int iy = 0; iy < size_F[1]; iy++)
	{
    
    int iy_glob = start_F[1] + iy;
    
    if(iy_glob < N/2)
    {
      ky[iy] = iy_glob * dk;
    }
    else
    {
      ky[iy] = (iy_glob-N) * dk;
		}
	}
  
  // kx
  for(int iz = 0; iz < size_F[2]; iz++)
	{
		int iz_glob = start_F[2] + iz;
		
    kz[iz] = iz_glob * dk;
	}
  
  // k2
  for(int ix = 0; ix<size_F[0]; ix++){
  for(int iy = 0; iy<size_F[1]; iy++){
  for(int iz = 0; iz<size_F[2]; iz++){

    int id = ix * size_F[1]*size_F[2] + iy * size_F[2] + iz;
    
    k2[id] = kx[ix]*kx[ix] + ky[iy]*ky[iy] + kz[iz]*kz[iz];
    
    k2[id] = std::max(k2[id], 1.e-12); 
    
  }}}
  
}

void CSpecDyn::setup_fields()
{
  switch(setup)
  {
    case 0:
      // u = 0, B = 0
      for(int id = 0; id < size_R_tot; id++)
      {    
        Vx_R[id] = 0.;
        Vy_R[id] = 0.;
        Vz_R[id] = 0.;
        
        Bx_R[id] = 0.;
        By_R[id] = 0.;
        Bz_R[id] = 0.;
      }
      break;
      
    case 1:
      // testing purposes
      for(int ix = 0; ix < size_R[0]; ix++){
      for(int iy = 0; iy < size_R[1]; iy++){
      for(int iz = 0; iz < size_R[2]; iz++){
      
        int id = ix * size_R[1]*size_R[2] + iy * size_R[2] + iz;
        
        double x_val = (start_R[0]+ix)*dx+XB;
        double y_val = (start_R[1]+iy)*dx+XB;
        double z_val = (start_R[2]+iz)*dx+XB;
        double r_val = sqrt(y_val*y_val + z_val*z_val);
        double p_val = atan2(z_val, y_val);
        double s_val = x_val;
        
        Vx_R[id] = 1.;
        Vy_R[id] = 2.;
        Vz_R[id] = 3.;
        
        Bx_R[id] = 0.;
        By_R[id] = 0.;
        Bz_R[id] = 0.;
    
      }}}
      break;
    
    default: 
      if(myRank==0){printf("No valid setup provided! setup = %d\n", setup);}
      MPI_Barrier(comm);
      MPI_Finalize();
      exit(EXIT_FAILURE);
  }
  
}

void CSpecDyn::execute()
{
  //~ for(int i = 0; i < 3; i++)
  //~ {
    print_vti();
    //~ time += dt;
  //~ }
}

void CSpecDyn::finalize()
{
  MPI_Finalize();
}

void CSpecDyn::print_vti()
{
  
  std::string file_name  = out_dir + "/step_" + std::to_string(vti_count) + ".vti";
  std::ofstream os;
  
  int offset = 0;
	int N_tot = N*N*N;
	int N_bytes_scalar  =   N_tot * sizeof(float);
	int N_bytes_vector  = 3*N_tot * sizeof(float);
  int bin_size_scalar = N_bytes_scalar + sizeof(uint64_t);// 2nd term is the size of the the leading integer announcing the numbers n the data chunk
  int bin_size_vector = N_bytes_vector + sizeof(uint64_t);
  
  // header
  if(myRank==0)
  {
    printf("Printing vti # %d!\n", vti_count);
    
    os.open(file_name.c_str(), std::ios::out);
    if(!os){
      std::cout << "Cannot write header to file '" << file_name << "'!\n";
    }
    
    // write header	
		int extend_l[3]  = {0, 0, 0};
		int extend_r[3]  = {N-1, N-1, N-1};
		double origin[3] = {XB,XB,XB};
    
    os << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;	
    os << "  <ImageData WholeExtent=\"" << extend_l[0] << " " << extend_r[0] << " " 
                                        << extend_l[1] << " " << extend_r[1] << " " 
                                        << extend_l[2] << " " << extend_r[2] 
				 << "\" Origin=\""  << origin[0]  << " " << origin[1]  << " " << origin[2] 
				 << "\" Spacing=\"" << dx << " " << dx << " " << dx
         << "\" Direction=\"0 0 1 0 1 0 1 0 0\">" << std::endl; // C -> FORTRAN order
    
    os << "      <FieldData>" << std::endl;
    os << "        <DataArray type=\"Float32\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">" << std::endl;
    os << "        "<< float(time) << std::endl;
    os << "        </DataArray>" << std::endl;
    os << "      </FieldData>" << std::endl;
        
		os << "    <Piece Extent=\"" << extend_l[0] << " " << extend_r[0] << " " 
                                 << extend_l[1] << " " << extend_r[1] << " " 
                                 << extend_l[2] << " " << extend_r[2] << "\">" << std::endl;
    
    os << "      <PointData Scalars=\"P\" Vectors=\"V\">" << std::endl;
    os << "        <DataArray type=\"Float32\" Name=\"Vx\" format=\"appended\" offset=\"" << offset << "\">" << std::endl;
    os << "        </DataArray>" << std::endl;
    offset += bin_size_scalar;
    os << "        <DataArray type=\"Float32\" Name=\"Vy\" format=\"appended\" offset=\"" << offset << "\">" << std::endl;
    os << "        </DataArray>" << std::endl;
    offset += bin_size_scalar;
    os << "        <DataArray type=\"Float32\" Name=\"Vz\" format=\"appended\" offset=\"" << offset << "\">" << std::endl;
    os << "        </DataArray>" << std::endl;
    offset += bin_size_scalar;
    os << "        <DataArray type=\"Float32\" Name=\"Bx\" format=\"appended\" offset=\"" << offset << "\">" << std::endl;
    os << "        </DataArray>" << std::endl;
    offset += bin_size_scalar;
    os << "        <DataArray type=\"Float32\" Name=\"By\" format=\"appended\" offset=\"" << offset << "\">" << std::endl;
    os << "        </DataArray>" << std::endl;
    offset += bin_size_scalar;
    os << "        <DataArray type=\"Float32\" Name=\"Bz\" format=\"appended\" offset=\"" << offset << "\">" << std::endl;
    os << "        </DataArray>" << std::endl;
    offset += bin_size_scalar;
    
    os << "      </PointData>" << std::endl;
    os << "      <CellData>" << std::endl;
    os << "      </CellData>" << std::endl;
    os << "    </Piece>" << std::endl;
    os << "  </ImageData>" << std::endl;
    os << "  <AppendedData encoding=\"raw\">" << std::endl;
    os << "   _";
                                
    os.close();
  
  }MPI_Barrier(comm);
  
  double start_time = MPI_Wtime();
  
  // print field components in parallel
  print_mpi_scalar(Vx_R, N_bytes_scalar, file_name.c_str());
  print_mpi_scalar(Vy_R, N_bytes_scalar, file_name.c_str()); 
  print_mpi_scalar(Vz_R, N_bytes_scalar, file_name.c_str()); 
  print_mpi_scalar(Bx_R, N_bytes_scalar, file_name.c_str());
  print_mpi_scalar(By_R, N_bytes_scalar, file_name.c_str()); 
  print_mpi_scalar(Bz_R, N_bytes_scalar, file_name.c_str()); 
  
  double print_time = MPI_Wtime() - start_time;
  if(myRank==0){printf("Print time = %f, pdims = [%d,%d], N = %d\n", print_time, pdims[0], pdims[1], N);}
  
  // footer
  if(myRank==0)
  {
    os.open(file_name.c_str(), std::ios::out | std::ios::app);
    if(!os){
      std::cout << "Cannot write footer to file '" << file_name << "'!\n";
      exit(3);
    }
		
		os << std::endl << "  </AppendedData>" << std::endl;
    os<< "</VTKFile>" << std::endl;
	
    os.close();
  }MPI_Barrier(comm);
  
  vti_count++;
}

void CSpecDyn::print_mpi_scalar(double* field, int& N_bytes_scalar, const char* file_name)
{
  if(myRank==0)
  {
    std::ofstream binary_os(file_name, std::ios::out | std::ios::app | std::ios::binary );
    binary_os.write(reinterpret_cast<const char*>(&N_bytes_scalar),sizeof(uint64_t)); // size of following binary package
    binary_os.close();
  }MPI_Barrier(comm);
  
  // open file
  MPI_File mpi_file;
  MPI_File_open(comm, file_name, MPI_MODE_APPEND|MPI_MODE_WRONLY, MPI_INFO_NULL, &mpi_file);
  
  // offset to end of file
  MPI_Offset mpi_eof;
  MPI_File_get_position(mpi_file, &mpi_eof);
  MPI_Barrier(comm);
  
  // data to float array
  for(int id = 0; id < size_R_tot; id++)
  {
    float_array[id] = float(field[id]);
  }
  
  // write data
  MPI_File_set_view(mpi_file, mpi_eof, MPI_FLOAT, vti_subarray, "native", MPI_INFO_NULL);
  MPI_File_write_all(mpi_file, float_array, size_R_tot, MPI_FLOAT, MPI_STATUS_IGNORE);
  
  // close file
  MPI_File_close(&mpi_file);  
}
