

//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "openclPort.h"
#include "opencl.h"

#include <fstream>
#include <sstream>
// #include <iostream>           // removed: avoid std::cerr / std::cout
#include <string>
#include <filesystem>  // C++17
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <R.h>                  // added: for Rprintf

namespace fs = std::filesystem;
using namespace openclPort;

// Load a single file like "nmath/bd0.cl"
namespace openclPort {

#ifdef USE_OPENCL
std::string load_kernel_source(const std::string& relative_path,
                               const std::string& package ) {
  // Retrieve full path via system.file()
  std::string path = Rcpp::as<std::string>(
    Rcpp::Function("system.file")("cl", relative_path,
                   Rcpp::Named("package") = package)
  );
  
  // Check for empty path returned by system.file (means file not found)
  if (path.empty()) {
    throw std::runtime_error("Kernel source not found via system.file: " + relative_path);
  }
  
  // Attempt to open the file
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("Failed to open kernel source: " + path);
  }
  
  // Read file contents
  std::ostringstream oss;
  oss << file.rdbuf();
  return oss.str();
}
#endif

/////////////////////////////

#ifdef USE_OPENCL
std::string load_kernel_library(const std::string& subdir, const std::string& package , bool verbose ) {
  std::string dir_path = Rcpp::as<std::string>(
    Rcpp::Function("system.file")("cl", subdir, Rcpp::Named("package") = package)
  );
  
  std::map<std::string, std::set<std::string>> provides_map;
  std::map<std::string, std::set<std::string>> depends_map;
  std::map<std::string, std::filesystem::path> file_map;
  
  if (verbose)  Rprintf("\n📂 Files found in '%s':\n", subdir.c_str());
  for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
    if (entry.path().extension() == ".cl") {
      std::string file_id = entry.path().stem().string();
      if (verbose) Rprintf(" - %s\n", file_id.c_str());
      
      std::ifstream infile(entry.path());
      std::string line;
      std::set<std::string> provides, depends;
      
      while (std::getline(infile, line)) {
        if (line.find("@provides") != std::string::npos) {
          std::stringstream ss(line.substr(line.find("@provides") + 9));
          std::string item;
          while (ss >> item) provides.insert(item);
        } else if (line.find("@depends") != std::string::npos) {
          std::stringstream ss(line.substr(line.find("@depends") + 9));
          std::string item;
          while (ss >> item) {
            // Remove only ‘,’ characters
            item.erase(std::remove(item.begin(), item.end(), ','), item.end());
            // item.erase(std::remove_if(item.begin(), item.end(), ::ispunct), item.end());
            depends.insert(item);
          }
        }
      }
      
      file_map[file_id] = entry.path();
      provides_map[file_id] = provides;
      depends_map[file_id] = depends;
    }
  }
  
  std::vector<std::string> sorted;
  std::set<std::string> sorted_set;
  std::set<std::string> unsorted_set;
  
  if (verbose)  Rprintf("\n📤 Files with no dependencies:\n");
  for (const auto& [file, _] : file_map) {
    if (depends_map[file].empty()) {
      sorted.push_back(file);
      sorted_set.insert(file);
      if (verbose) Rprintf(" + %s\n", file.c_str());
    } else {
      unsorted_set.insert(file);
    }
  }
  
  if (verbose)  Rprintf("\n🧪 Unsorted files:\n");
  for (const auto& file : unsorted_set) {
    if (verbose) Rprintf(" - %s\n", file.c_str());
  }
  
  int pass_count = 0;
  while (!unsorted_set.empty()) {
    ++pass_count;
    if (verbose) Rprintf("\n🔁 While Loop Pass #%d — Remaining unsorted: %d\n", pass_count, (int)unsorted_set.size());
    
    std::vector<std::string> newly_sorted;
    bool progress_made = false;
    int file_counter = 0;
    
    for (const std::string& file : unsorted_set) {
      ++file_counter;
      if (verbose) Rprintf("   🔍 File #%d: %s\n", file_counter, file.c_str());
      
      const auto& deps = depends_map[file];
      int depends_counter = static_cast<int>(deps.size());
      if (verbose) Rprintf("      📦 Dependency Count: %d\n", depends_counter);
      
      int found_counter = 0;
      int dep_index = 0;
      for (const std::string& dep : deps) {
        ++dep_index;
        if (verbose) Rprintf("         🔎 Checking classified #%d: %s\n", dep_index, dep.c_str());
        
        auto it = sorted_set.find(dep);
        if (it != sorted_set.end()) {
          if (verbose) Rprintf("            ➤ Found in sorted? ✅ Yes\n");
          ++found_counter;
        } else {
          if (verbose) Rprintf("            ➤ Found in sorted? ❌ No\n");
        }
      }
      
      if (verbose) Rprintf("      🔍 Found count: %d\n", found_counter);
      if (found_counter == depends_counter) {
        sorted.push_back(file);
        sorted_set.insert(file);
        newly_sorted.push_back(file);
        progress_made = true;
        if (verbose) Rprintf(" ✅ Promoted to Sorted: %s\n", file.c_str());
      }
    }
    
    for (const std::string& file : newly_sorted) {
      unsorted_set.erase(file);
    }
    
    if (!progress_made) {
      if (verbose) {
        Rprintf("\n❌ No files promoted on pass #%d; possible circular or missing dependencies:\n", pass_count);
        for (const std::string& file : unsorted_set) {
          Rprintf(" - %s\n", file.c_str());
        }
      }
      throw std::runtime_error("Dependency sort failed: unresolved dependencies remain.");
    }
  }
  
  if (verbose)  Rprintf("\n🔗 Final Sorted Load Order:\n");
  for (const auto& file : sorted) {
    if (verbose) Rprintf(" - %s\n", file.c_str());
  }
  
  std::string combined_source;
  for (const auto& file : sorted) {
    std::string rel_path = subdir + "/" + file + ".cl";
    combined_source += load_kernel_source(rel_path, package) + "\n";
  }
  
  return combined_source;
}

namespace {

std::vector<std::string> parse_cl_tag(
    const std::vector<std::string>& lines,
    const std::string& tag)
{
  std::vector<std::string> result;
  std::string pattern = "@" + tag;
  for (const auto& line : lines) {
    auto pos = line.find(pattern);
    if (pos == std::string::npos) continue;
    auto colon = line.find(':', pos + pattern.size());
    if (colon == std::string::npos) continue;
    std::istringstream ss(line.substr(colon + 1));
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      tok.erase(0, tok.find_first_not_of(" \t\r\n"));
      auto last = tok.find_last_not_of(" \t\r\n");
      if (last != std::string::npos) tok.erase(last + 1);
      if (!tok.empty()) result.push_back(tok);
    }
  }
  return result;
}

struct KernelDepIndex {
  std::vector<std::string> stems_ordered;
  std::unordered_map<std::string, std::vector<std::string>> all_depends;
};

KernelDepIndex read_tsv_index(const std::string& tsv_path)
{
  KernelDepIndex idx;
  std::ifstream f(tsv_path);
  if (!f.is_open()) {
    throw std::runtime_error(
        "kernel_dependency_index.tsv not found: " + tsv_path);
  }
  std::string line;
  bool header = true;
  while (std::getline(f, line)) {
    if (header) { header = false; continue; }
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;

    auto tab = line.find('\t');
    std::string stem = (tab == std::string::npos) ? line : line.substr(0, tab);
    if (stem.empty()) continue;
    idx.stems_ordered.push_back(stem);

    std::vector<std::string> deps;
    if (tab != std::string::npos && tab + 1 < line.size()) {
      std::istringstream ss(line.substr(tab + 1));
      std::string tok;
      while (std::getline(ss, tok, ',')) {
        tok.erase(0, tok.find_first_not_of(" \t\r\n"));
        auto last = tok.find_last_not_of(" \t\r\n");
        if (last != std::string::npos) tok.erase(last + 1);
        if (!tok.empty()) deps.push_back(tok);
      }
    }
    idx.all_depends[stem] = std::move(deps);
  }
  return idx;
}

std::string load_library_for_kernel_cross_package(
    const std::string& kernel_relative_path,
    const std::string& kernel_package,
    const std::string& library_subdir,
    const std::string& library_package,
    const std::string& depends_tag)
{
  std::string kernel_path = Rcpp::as<std::string>(
      Rcpp::Function("system.file")(
          "cl", kernel_relative_path,
          Rcpp::Named("package") = kernel_package));
  if (kernel_path.empty()) {
    throw std::runtime_error(
        "Kernel file not found via system.file: " + kernel_relative_path +
        " (package=" + kernel_package + ")");
  }

  std::string lib_dir = Rcpp::as<std::string>(
      Rcpp::Function("system.file")(
          "cl", library_subdir,
          Rcpp::Named("package") = library_package));
  if (lib_dir.empty()) {
    throw std::runtime_error(
        "Library directory not found via system.file: " + library_subdir +
        " (package=" + library_package + ")");
  }

  std::ifstream kf(kernel_path);
  if (!kf.is_open()) {
    throw std::runtime_error("Cannot open kernel file: " + kernel_path);
  }
  std::vector<std::string> klines;
  {
    std::string kl;
    while (std::getline(kf, kl)) klines.push_back(kl);
  }
  kf.close();

  std::vector<std::string> needed_stems = parse_cl_tag(klines, depends_tag);
  if (needed_stems.empty()) {
    return "";
  }

  std::string tsv_path = lib_dir + "/kernel_dependency_index.tsv";
  KernelDepIndex idx = read_tsv_index(tsv_path);

  std::unordered_set<std::string> needed_set(needed_stems.begin(), needed_stems.end());
  std::vector<std::string> to_load;
  to_load.reserve(needed_set.size());
  for (const auto& stem : idx.stems_ordered) {
    if (needed_set.count(stem)) to_load.push_back(stem);
  }

  std::string combined;
  for (const auto& stem : to_load) {
    std::string cl_path = lib_dir + "/" + stem + ".cl";
    std::ifstream cf(cl_path);
    if (!cf.is_open()) {
      throw std::runtime_error(
          "Library file not found for stem '" + stem + "': " + cl_path);
    }
    std::ostringstream oss;
    oss << cf.rdbuf();
    combined += oss.str() + "\n\n";
  }

  return combined;
}

std::string load_library_for_kernel(
    const std::string& kernel_relative_path,
    const std::string& library_subdir,
    const std::string& package,
    const std::string& depends_tag)
{
  return load_library_for_kernel_cross_package(
      kernel_relative_path,
      package,
      library_subdir,
      package,
      depends_tag);
}

std::string resolve_kernel_path(
    const std::string& family,
    const std::string& link)
{
  if (family == "binomial" || family == "quasibinomial") {
    if (link == "logit") {
      return "src/f2_f3_binomial_logit.cl";
    }
    if (link == "probit") {
      return "src/f2_f3_binomial_probit.cl";
    }
    if (link == "cloglog") {
      return "src/f2_f3_binomial_cloglog.cl";
    }
    throw std::runtime_error(
        "Unsupported link function for binomial family: " + link);
  }
  if (family == "poisson" || family == "quasipoisson") {
    return "src/f2_f3_poisson.cl";
  }
  if (family == "Gamma") {
    return "src/f2_f3_gamma.cl";
  }
  if (family == "gaussian") {
    return "src/f2_f3_gaussian.cl";
  }
  throw std::runtime_error("Unsupported family: " + family);
}

} // namespace

#endif // USE_OPENCL

} // namespace openclPort



#ifdef USE_OPENCL
namespace glmbayes {
namespace opencl {

std::string load_likelihood_subgradient_program(
    const std::string& family,
    const std::string& link,
    const std::string& package)
{
  const std::string kernel_file = resolve_kernel_path(family, link);

  std::string opencl_source = openclPort::load_kernel_source("OPENCL.cl", package);
  std::string libr_shims_source =
      openclPort::load_kernel_library("libR_shims", package, false);
  std::string r_ext_types_source =
      openclPort::load_kernel_library("R_ext_types", package, false);
  std::string r_shims_source =
      openclPort::load_kernel_library("R_shims", package, false);
  std::string r_ext_runtime_source =
      openclPort::load_kernel_library("R_ext_runtime", package, false);
  std::string r_ext_internals_source =
      openclPort::load_kernel_library("R_ext_internals", package, false);
  std::string system_source =
      openclPort::load_kernel_library("System", package, false);
  std::string nmath_source = load_library_for_kernel(
      kernel_file, "nmath", package, "all_depends_nmath");
  std::string ksrc = openclPort::load_kernel_source(kernel_file, package);

  return opencl_source +
    "\n" + libr_shims_source +
    "\n" + r_ext_types_source +
    "\n" + r_shims_source +
    "\n" + r_ext_runtime_source +
    "\n" + r_ext_internals_source +
    "\n" + system_source +
    "\n" + nmath_source +
    "\n" + ksrc;
}

std::string load_likelihood_subgradient_program_v2(
    const std::string& family,
    const std::string& link,
    const std::string& app_package,
    const std::string& nmath_package)
{
  const std::string kernel_file = resolve_kernel_path(family, link);

  std::string opencl_source =
      openclPort::load_kernel_source("OPENCL.cl", nmath_package);
  std::string libr_shims_source =
      openclPort::load_kernel_library("libR_shims", nmath_package, false);
  std::string r_ext_types_source =
      openclPort::load_kernel_library("R_ext_types", nmath_package, false);
  std::string r_shims_source =
      openclPort::load_kernel_library("R_shims", nmath_package, false);
  std::string r_ext_runtime_source =
      openclPort::load_kernel_library("R_ext_runtime", nmath_package, false);
  std::string r_ext_internals_source =
      openclPort::load_kernel_library("R_ext_internals", nmath_package, false);
  std::string system_source =
      openclPort::load_kernel_library("System", nmath_package, false);
  std::string nmath_source = load_library_for_kernel_cross_package(
      kernel_file,
      app_package,
      "nmath",
      nmath_package,
      "all_depends_nmath");
  std::string ksrc =
      openclPort::load_kernel_source(kernel_file, app_package);

  return opencl_source +
    "\n" + libr_shims_source +
    "\n" + r_ext_types_source +
    "\n" + r_shims_source +
    "\n" + r_ext_runtime_source +
    "\n" + r_ext_internals_source +
    "\n" + system_source +
    "\n" + nmath_source +
    "\n" + ksrc;
}

} // namespace opencl
} // namespace glmbayes
#endif // USE_OPENCL

