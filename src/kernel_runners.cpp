

#ifdef USE_OPENCL
#include "kernel_loader.h"

#ifdef USE_DIRECT_CLH
// we passed “-I…/include/CL -DUSE_DIRECT_CLH”
#include <CL/cl.h>
#else
// normal case on Linux/macOS/Windows
#include <CL/cl.h>
#endif
#endif


//#include <Rcpp.h>
# include <RcppArmadillo.h>
#include <vector>
#include <string>


#ifdef USE_OPENCL

void f2_binomial_logit_prep_kernel_runner(
    const std::string& kernel_source,
    const char*        kernel_name,
    int                l1,
    int                l2,
    int                m1,
    const std::vector<double>& X_flat,
    const std::vector<double>& B_flat,
    const std::vector<double>& mu_flat,
    const std::vector<double>& P_flat,
    const std::vector<double>& alpha_flat,
    std::vector<double>&       qf_flat,   // output, length = m1
    std::vector<double>&       xb_flat,   // output, length = l1*m1
    int progbar 
) {
  if ((int)X_flat.size() != l1*l2 ||
      (int)B_flat.size() != l2*m1 ||
      (int)mu_flat.size() != l2 ||
      (int)P_flat.size() != l2*l2 ||
      (int)alpha_flat.size() != l1) {
    throw std::runtime_error("Input flat-vector sizes mismatch dimensions.");
  }
  
  // Initialize outputs
  qf_flat.assign(m1, 0.0);
  xb_flat.assign((size_t)l1*m1, 0.0);
  
  cl_int status;
  // 🔍 Platform & Device
  cl_platform_id platform;
  cl_device_id   device;
  status = clGetPlatformIDs(1, &platform, nullptr);
  status |= clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
  
  // 🌐 Context & Queue
  cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  cl_queue_properties props[] = {0};
  cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, props, &status);
  
  // 📦 Program & Kernel
  const char* src_ptr = kernel_source.c_str();
  size_t src_len = kernel_source.size();
  cl_program program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  
  // 🧮 Device Buffers
  cl_mem bufX = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                               sizeof(double)*X_flat.size(), (void*)X_flat.data(), &status);
  cl_mem bufB = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                               sizeof(double)*B_flat.size(), (void*)B_flat.data(), &status);
  cl_mem bufMu= clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                               sizeof(double)*mu_flat.size(), (void*)mu_flat.data(), &status);
  cl_mem bufP = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                               sizeof(double)*P_flat.size(), (void*)P_flat.data(), &status);
  cl_mem bufA = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                               sizeof(double)*alpha_flat.size(), (void*)alpha_flat.data(), &status);
  
  cl_mem bufQF= clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                               sizeof(double)*qf_flat.size(), nullptr, &status);
  cl_mem bufXB= clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                               sizeof(double)*xb_flat.size(), nullptr, &status);
  
  // ⏫ Transfer Data to Device (already via COPY_HOST_PTR for inputs)
  
  // 🗳️ Set Kernel Args
  int arg = 0;
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufX);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufB);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufMu);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufP);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufA);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufQF);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufXB);
  clSetKernelArg(kernel, arg++, sizeof(int), &l1);
  clSetKernelArg(kernel, arg++, sizeof(int), &l2);
  clSetKernelArg(kernel, arg++, sizeof(int), &m1);
  
  // 🚀 Launch Parallel Kernel
  size_t global = m1;
  
//  size_t global = m1;
  
  status = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
  
  // 📥 Retrieve Output
  status = clEnqueueReadBuffer(queue, bufQF, CL_TRUE, 0,
                               sizeof(double)*qf_flat.size(), qf_flat.data(),
                               0, nullptr, nullptr);
  status = clEnqueueReadBuffer(queue, bufXB, CL_TRUE, 0,
                               sizeof(double)*xb_flat.size(), xb_flat.data(),
                               0, nullptr, nullptr);
  
  // 🧹 Cleanup
  clReleaseMemObject(bufX);
  clReleaseMemObject(bufB);
  clReleaseMemObject(bufMu);
  clReleaseMemObject(bufP);
  clReleaseMemObject(bufA);
  clReleaseMemObject(bufQF);
  clReleaseMemObject(bufXB);
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(queue);
  clReleaseContext(context);
}
#endif



#ifdef USE_OPENCL

void f2_binomial_logit_prep_grad_kernel_runner(
    const std::string&        kernel_source,  // your .cl contents
    const char*               kernel_name,    // "f2_binomial_logit_prep_grad"
    int                       l1,             // nobs
    int                       l2,             // ncoef
    int                       m1,             // ngrids
    const std::vector<double>& X_flat,        // length = l1*l2
    const std::vector<double>& B_flat,        // length = m1*l2
    const std::vector<double>& mu_flat,       // length = l2
    const std::vector<double>& P_flat,        // length = l2*l2
    const std::vector<double>& alpha_flat,    // length = l1
    const std::vector<double>& y_flat,        // length = l1
    const std::vector<double>& wt_flat,       // length = l1
    std::vector<double>&       qf_flat,       // OUT: length = m1
    std::vector<double>&       xb_flat,       // OUT: length = m1*l1
    std::vector<double>&       grad_flat,     // OUT: length = m1*l2
    int                       progbar 
) {
  // 0) Sanity‐check sizes
  if ((int)X_flat.size()   != l1*l2 ||
      (int)B_flat.size()   != m1*l2 ||
      (int)mu_flat.size()  != l2 ||
      (int)P_flat.size()   != l2*l2 ||
      (int)alpha_flat.size()!= l1 ||
      (int)y_flat.size()   != l1 ||
      (int)wt_flat.size()  != l1) {
    throw std::runtime_error("Input flat‐vector sizes mismatch dimensions.");
  }
  
  // 1) Initialize output buffers
  qf_flat  .assign(m1,           0.0);
  xb_flat  .assign((size_t)l1*m1, 0.0);
  grad_flat.assign((size_t)l2*m1, 0.0);
  
  cl_int status;
  
  // 2) Platform & Device
  cl_platform_id platform;
  cl_device_id   device;
  status  = clGetPlatformIDs(1, &platform, nullptr);
  status |= clGetDeviceIDs  (platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);
  
  // 3) Context & Queue
  cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  cl_queue_properties props[] = {0};
  cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, props, &status);
  
  // 4) Program & Kernel
  const char* src_ptr = kernel_source.c_str();
  size_t      src_len = kernel_source.size();
  cl_program  program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  status |= clBuildProgram   (program, 0, nullptr, nullptr, nullptr, nullptr);
  
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  
  // 5) Device Buffers
  cl_mem bufX    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*X_flat.size(),   (void*)X_flat.data(),   &status);
  cl_mem bufB    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*B_flat.size(),   (void*)B_flat.data(),   &status);
  cl_mem bufMu   = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*mu_flat.size(),  (void*)mu_flat.data(),  &status);
  cl_mem bufP    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*P_flat.size(),   (void*)P_flat.data(),   &status);
  cl_mem bufA    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*alpha_flat.size(), (void*)alpha_flat.data(), &status);
  cl_mem bufY    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*y_flat.size(),   (void*)y_flat.data(),   &status);
  cl_mem bufW    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*wt_flat.size(),  (void*)wt_flat.data(),  &status);
  
  cl_mem bufQF   = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                  sizeof(double)*qf_flat.size(),   nullptr, &status);
  cl_mem bufXB   = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                  sizeof(double)*xb_flat.size(),   nullptr, &status);
  cl_mem bufGrad = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                  sizeof(double)*grad_flat.size(), nullptr, &status);
  
  // 6) Set Kernel Args (must match the .cl signature exactly)
  int arg = 0;
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufX);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufB);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufMu);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufP);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufA);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufY);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufW);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufQF);
//  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufXB);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufGrad);
  clSetKernelArg(kernel, arg++, sizeof(int),    &l1);
  clSetKernelArg(kernel, arg++, sizeof(int),    &l2);
  clSetKernelArg(kernel, arg++, sizeof(int),    &m1);
  
  // 7) Launch
  size_t global = (size_t)m1;
  status = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
  
  // 8) Read back outputs
  status = clEnqueueReadBuffer(queue, bufQF,   CL_TRUE, 0,
                               sizeof(double)*qf_flat.size(),   qf_flat.data(),
                               0, nullptr, nullptr);
  
  status = clEnqueueReadBuffer(queue, bufXB,   CL_TRUE, 0,
                               sizeof(double)*xb_flat.size(),   xb_flat.data(),
                               0, nullptr, nullptr);
  
  status = clEnqueueReadBuffer(queue, bufGrad, CL_TRUE, 0,
                               sizeof(double)*grad_flat.size(), grad_flat.data(),
                               0, nullptr, nullptr);
  
  
  // 8a) Sanity-check: error out if both outputs are all zeros
                               {
                                 auto all_zero = [](auto& vec){
                                   return std::all_of(vec.begin(), vec.end(),
                                                      [](double x){ return x == 0.0; });
                                 };
                                 
                                 bool qf_is_zero   = all_zero(qf_flat);
                                 bool grad_is_zero = all_zero(grad_flat);
                                 
                                 if (qf_is_zero || grad_is_zero) {
                                   std::ostringstream msg;
                                   msg << "OpenCL kernel returned "
                                       << (qf_is_zero   ? "qf_flat all zeros "   : "")
                                       << (grad_is_zero ? "grad_flat all zeros." : "");
                                   throw std::runtime_error(msg.str());
                                 }
                               }
  
  // --- Begin modified cleanup section ---
  // 9a) Drain any pending commands
  clFlush(queue);
  clFinish(queue);
  
  // 9b) Release buffers (inverse creation order)
  clReleaseMemObject(bufGrad);
  clReleaseMemObject(bufXB);
  clReleaseMemObject(bufQF);
  clReleaseMemObject(bufW);
  clReleaseMemObject(bufY);
  clReleaseMemObject(bufA);
  clReleaseMemObject(bufP);
  clReleaseMemObject(bufMu);
  clReleaseMemObject(bufB);
  clReleaseMemObject(bufX);
  
  // 9c) Release kernel, program, queue, context
  clReleaseKernel       (kernel);
  clReleaseProgram      (program);
  clReleaseCommandQueue (queue);
  clReleaseContext      (context);
  // --- End modified cleanup section ---
  
  }
#endif



#ifdef USE_OPENCL

void f2_f3_kernel_runner(
    const std::string&        kernel_source,
    const char*               kernel_name,
    int                       l1,
    int                       l2,
    int                       m1,
    const std::vector<double>& X_flat,
    const std::vector<double>& B_flat,
    const std::vector<double>& mu_flat,
    const std::vector<double>& P_flat,
    const std::vector<double>& alpha_flat,
    const std::vector<double>& y_flat,
    const std::vector<double>& wt_flat,
    std::vector<double>&       qf_flat,
    std::vector<double>&       grad_flat,
    int                        progbar
) {
  // 0) Sanity-check sizes
  if ((int)X_flat.size()    != l1*l2 ||
      (int)B_flat.size()    != m1*l2 ||
      (int)mu_flat.size()   != l2    ||
      (int)P_flat.size()    != l2*l2 ||
      (int)alpha_flat.size()!= l1    ||
      (int)y_flat.size()    != l1    ||
      (int)wt_flat.size()   != l1) {
    throw std::runtime_error("Input flat-vector sizes mismatch dimensions.");
  }
  
  // 1) Initialize output buffers
  qf_flat.assign(m1, 0.0);
  grad_flat.assign((size_t)l2*m1, 0.0);
  
  cl_int status = 0;
  
  // 2) Platform & Device
  Rcpp::Rcout << "[runner] P0: before clGetPlatformIDs\n";
  
  cl_platform_id platform = nullptr;
  cl_device_id   device   = nullptr;
  
  // ---- Platform ----
  status = clGetPlatformIDs(1, &platform, nullptr);

  // -1001 = CL_PLATFORM_NOT_FOUND_KHR (not always defined in headers)
  if (status == -1001) {
    throw std::runtime_error(
        "OpenCL error: no OpenCL platforms found (clGetPlatformIDs returned -1001). "
        "Your system does not expose an OpenCL platform."
    );
  }
  if (status != CL_SUCCESS) {
    std::ostringstream msg;
    msg << "OpenCL error: clGetPlatformIDs failed with status " << status << ".";
    throw std::runtime_error(msg.str());
  }
  
  // ---- Device ----
  status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_DEFAULT, 1, &device, nullptr);

  // -9 = CL_DEVICE_NOT_FOUND
  if (status == -9) {
    throw std::runtime_error(
        "OpenCL error: no suitable OpenCL GPU devices found "
        "(clGetDeviceIDs returned -9)."
    );
  }
  if (status != CL_SUCCESS) {
    std::ostringstream msg;
    msg << "OpenCL error: clGetDeviceIDs failed with status " << status << ".";
    throw std::runtime_error(msg.str());
  }
  
  // 3) Context & Queue
  Rcpp::Rcout << "[runner] P3: before clCreateContext\n";
  cl_context context = clCreateContext(nullptr, 1, &device, nullptr, nullptr, &status);
  Rcpp::Rcout << "[runner] P4: after clCreateContext, status=" << status << "\n";
  
  cl_queue_properties props[] = {0};
  Rcpp::Rcout << "[runner] P5: before clCreateCommandQueueWithProperties\n";
  cl_command_queue queue = clCreateCommandQueueWithProperties(context, device, props, &status);
  Rcpp::Rcout << "[runner] P6: after clCreateCommandQueueWithProperties, status=" << status << "\n";
  
  // 4) Program & Kernel
  const char* src_ptr = kernel_source.c_str();
  size_t      src_len = kernel_source.size();
  Rcpp::Rcout << "[runner] P7: before clCreateProgramWithSource, src_len=" << src_len << "\n";
  cl_program  program = clCreateProgramWithSource(context, 1, &src_ptr, &src_len, &status);
  Rcpp::Rcout << "[runner] P8: after clCreateProgramWithSource, status=" << status << "\n";
  
  Rcpp::Rcout << "[runner] P9: before clBuildProgram\n";
  status |= clBuildProgram(program, 0, nullptr, nullptr, nullptr, nullptr);
  Rcpp::Rcout << "[runner] P10: after clBuildProgram, status=" << status << "\n";
  
  Rcpp::Rcout << "[runner] P11: before clCreateKernel\n";
  cl_kernel kernel = clCreateKernel(program, kernel_name, &status);
  Rcpp::Rcout << "[runner] P12: after clCreateKernel, status=" << status << "\n";
  
  // 5) Device Buffers
  Rcpp::Rcout << "[runner] A: before buffer creation\n";
  
  cl_mem bufX    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*X_flat.size(),   (void*)X_flat.data(),   &status);
  cl_mem bufB    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*B_flat.size(),   (void*)B_flat.data(),   &status);
  cl_mem bufMu   = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*mu_flat.size(),  (void*)mu_flat.data(),  &status);
  cl_mem bufP    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*P_flat.size(),   (void*)P_flat.data(),   &status);
  cl_mem bufA    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*alpha_flat.size(), (void*)alpha_flat.data(), &status);
  cl_mem bufY    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*y_flat.size(),   (void*)y_flat.data(),   &status);
  cl_mem bufW    = clCreateBuffer(context, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR,
                                  sizeof(double)*wt_flat.size(),  (void*)wt_flat.data(),  &status);
  
  cl_mem bufQF   = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                  sizeof(double)*qf_flat.size(),   nullptr, &status);
  cl_mem bufGrad = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                  sizeof(double)*grad_flat.size(), nullptr, &status);
  
  Rcpp::Rcout << "[runner] B: after buffer creation\n";
  
  // 6) Set Kernel Args
  int arg = 0;
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufX);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufB);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufMu);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufP);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufA);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufY);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufW);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufQF);
  clSetKernelArg(kernel, arg++, sizeof(cl_mem), &bufGrad);
  clSetKernelArg(kernel, arg++, sizeof(int),    &l1);
  clSetKernelArg(kernel, arg++, sizeof(int),    &l2);
  clSetKernelArg(kernel, arg++, sizeof(int),    &m1);
  
  // 7) Launch
  size_t global = (size_t)m1;
  Rcpp::Rcout << "[runner] C: before enqueue\n";
  status = clEnqueueNDRangeKernel(queue, kernel, 1, nullptr, &global, nullptr, 0, nullptr, nullptr);
  Rcpp::Rcout << "[runner] D: after enqueue\n";
  
  // 8) Read back outputs
  Rcpp::Rcout << "[runner] E: before read qf\n";
  status = clEnqueueReadBuffer(queue, bufQF,   CL_TRUE, 0,
                               sizeof(double)*qf_flat.size(),   qf_flat.data(),
                               0, nullptr, nullptr);
  Rcpp::Rcout << "[runner] F: after read qf\n";
  
  Rcpp::Rcout << "[runner] G: before read grad\n";
  status = clEnqueueReadBuffer(queue, bufGrad, CL_TRUE, 0,
                               sizeof(double)*grad_flat.size(), grad_flat.data(),
                               0, nullptr, nullptr);
  Rcpp::Rcout << "[runner] H: after read grad\n";
  
  // 8a) Sanity-check: error out if both outputs are all zeros
  {
    auto all_zero = [](auto& vec){
      return std::all_of(vec.begin(), vec.end(),
                         [](double x){ return x == 0.0; });
    };
    
    bool qf_is_zero   = all_zero(qf_flat);
    bool grad_is_zero = all_zero(grad_flat);
    
    if (qf_is_zero || grad_is_zero) {
      std::ostringstream msg;
      msg << "OpenCL kernel returned "
          << (qf_is_zero   ? "qf_flat all zeros "   : "")
          << (grad_is_zero ? "grad_flat all zeros." : "");
      throw std::runtime_error(msg.str());
    }
  }
  
  clFlush(queue);
  clFinish(queue);
  
  clReleaseMemObject(bufGrad);
  clReleaseMemObject(bufQF);
  clReleaseMemObject(bufW);
  clReleaseMemObject(bufY);
  clReleaseMemObject(bufA);
  clReleaseMemObject(bufP);
  clReleaseMemObject(bufMu);
  clReleaseMemObject(bufB);
  clReleaseMemObject(bufX);
  
  clReleaseKernel       (kernel);
  clReleaseProgram      (program);
  clReleaseCommandQueue (queue);
  clReleaseContext      (context);
}

#endif

#ifdef USE_OPENCL

int detect_num_gpus_internal() {
  cl_uint num_platforms = 0;
  cl_int status = clGetPlatformIDs(0, nullptr, &num_platforms);
  if (status != CL_SUCCESS || num_platforms == 0) return 0;
  
  std::vector<cl_platform_id> platforms(num_platforms);
  clGetPlatformIDs(num_platforms, platforms.data(), nullptr);
  
  int total_compute_units = 0;
  for (auto& platform : platforms) {
    cl_uint num_devices = 0;
    status = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, nullptr, &num_devices);
    if (status != CL_SUCCESS || num_devices == 0) continue;
    
    std::vector<cl_device_id> devices(num_devices);
    clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, num_devices, devices.data(), nullptr);
    
    for (auto& device : devices) {
      cl_uint units = 0;
      status = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(units), &units, nullptr);
      if (status == CL_SUCCESS) total_compute_units += units;
    }
  }
  return total_compute_units;
}

#endif