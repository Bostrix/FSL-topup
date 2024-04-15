# TOPUP (Tool for Estimating and Correcting Susceptibility-Induced Distortions) - Compilation and Usage Guide

# Introduction
Welcome to TOPUP, an essential tool for estimating and correcting susceptibility-induced distortions in diffusion imaging. This comprehensive guide will assist you in compiling the TOPUP tool and effectively utilizing its features for distortion correction in your imaging data.


For more information about TOPUP and related tools, visit the FMRIB Software Library (FSL) website: [FSL Git Repository](https://git.fmrib.ox.ac.uk/fsl)
You can also find additional resources and documentation on TOPUP on the FSL wiki page: [TOPUP Documentation](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup).

## Clone the Repository
Begin by cloning the project repository from GitHub onto your local machine. You can do this by running the following command in your terminal or command prompt:
```bash
git clone https://github.com/Bostrix/FSL-topup.git
```
This command will create a local copy of the project in a directory named "FSL-topup".

## Navigate to Project Directory
Change your current directory to the newly cloned project directory using the following command:
```bash
cd FSL-topup
```
## Installing Dependencies
To install the necessary dependencies for compiling and building the project, follow these steps:
```bash
sudo apt-get update
sudo apt install g++
sudo apt install make
sudo apt-get install libblas-dev libblas3
sudo apt-get install liblapack-dev liblapack3
sudo apt-get install zlib1g zlib1g-dev
sudo apt-get install libboost-all-dev
```
After completing these steps, you should have all the necessary dependencies installed on your system to use TOPUP.

## Compilation
To compile TOPUP, follow these steps:

- Ensure correct path in Makefile:
 After installing the necessary tools, verify correct path in the makefile to include additional LDFLAGS for the required libraries. For instance, if utilizing the warpfns library, basisfield library, meshclass library, miscmaths library,utils and znzlib, ensure that the correct path is present in the makefile.
Make sure `$(ADDED_LDFLAGS)` is included in the compile step of the makefile.

- Confirm that the file `point_list.h` within the warpfns library accurately includes the path to `armawrap/newmat.h`.
  
- Verify the accurate paths in meshclass's Makefile:
verify the correct path in the makefile of meshclass to include Linker flags for the required libraries. For instance, if utilizing the newimage,miscmaths,NewNifti,cprob,znzlib,utils libraries as LDFLAGS variable in meshclass makefile,ensure that the correct path is present in the makefile.

- Compile Source Code:
Execute the appropriate compile command to build the TOPUP tool. For example:
```bash
make clean
make
```
This command will compile the source code and generate the executable file for topup.

- Resolving Shared Library Errors:
When running an executable on Linux, you may encounter errors related to missing shared libraries.This typically manifests as messages like:
```bash
./topup: error while loading shared libraries: libexample.so: cannot open shared object file:No such file or directory
./applytopup: error while loading shared libraries: libexample.so: cannot open shared object file:No such file or directory
```
To resolve these errors,Pay attention to the names of the missing libraries mentioned in the error message.Locate the missing libraries on your system. If they are not present, you may need to install the corresponding packages.If the libraries are installed in custom directories, you need to specify those directories using the `LD_LIBRARY_PATH` environment variable. For example:
```bash
export LD_LIBRARY_PATH=/path/to/custom/libraries:$LD_LIBRARY_PATH
```
Replace `/path/to/custom/libraries` with the actual path to the directory containing the missing libraries.Once the `LD_LIBRARY_PATH` is set, attempt to run the executable again.If you encounter additional missing library errors, repeat steps until all dependencies are resolved.

- Resolving "The environment variable FSLOUTPUTTYPE is not defined" errors:
If you encounter an error related to the FSLOUTPUTTYPE environment variable not being set.Setting it to `NIFTI_GZ` is a correct solution, as it defines the output format for FSL tools to NIFTI compressed with gzip.Here's how you can resolve:
```bash
export FSLOUTPUTTYPE=NIFTI_GZ
```
By running this command, you've set the `FSLOUTPUTTYPE` environment variable to `NIFTI_GZ`,which should resolve the error you encountered.

## Usage
Once TOPUP is successfully compiled, you can use it for distortion correction in diffusion imaging. Follow these steps to utilize TOPUP and applytopup effectively:
- Running TOPUP
  Ensure you have prepared Input images with different acquisition parameters,Text file containing acquisition parameters and Configuration file specifying TOPUP parameters (optional but recommended).
  Execute the ./topup command with the required options and input files. Here is the basic syntax:
```bash
  topup --imain=<input_images.nii> --datain=<acquisition_parameters.txt> --config=<config_file.cnf> --out=<output_prefix>
```
  Replace `<input_images.nii>`, `<acquisition_parameters.txt>`, `<config_file.cnf>`, and `<output_prefix>` with the paths to your input images, acquisition parameters file, configuration file, and output prefix, respectively.
Customize the behavior of TOPUP by providing additional options as needed. Refer to the usage guide provided in the documentation for a list of available options and their descriptions.

- Applying TOPUP Corrections
Ensure you have prepared Input images to be corrected (e.g., b=0 and diffusion weighted images) and Output directory for corrected images.Then,Run applytopup with the following command syntax:
```bash
applytopup -i=<input_images.nii> -a=<acquisition_parameters.txt> -x=<index.txt> -t=<topup_results> -o=<output_directory>
```
Replace `<input_images.nii>`, `<acquisition_parameters.txt>`, `<index.txt>`, `<topup_results>`, and `<output_directory>` with the paths to your input images, acquisition parameters file, index file, TOPUP results, and output directory, respectively.

- Interpreting Output
After running TOPUP and applytopup, you will obtain various output files that provide valuable information about the distortion correction process.Customize the behavior of TOPUP by providing additional options as needed. Refer to the usage guide provided in the documentation for a list of available options and their descriptions.

# Conclusion
 You have now successfully compiled TOPUP and gained insights into its usage for distortion correction in diffusion imaging. By following the steps outlined in this guide, you can effectively utilize TOPUP to correct susceptibility-induced distortions in your imaging data.If you encounter any issues or have further questions, refer to the provided documentation or seek assistance from the project maintainers.
