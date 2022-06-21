#' Title
#'
#' @param count_table
#' @param cell_cluster_conversion
#' @param relative_prop
#' @param sample_new_levels
#' @param use_average_cluster_profiles
#' @param rate
#' @param cluster_size_max
#' @param cluster_size_min
#' @param sampling_type
#' @param nCV
#' @param simulation_parameter
#' @param simulation_model
#'
#' @return
#' @export
#'
#' @examples
sc2spatial=function(count_table, cell_cluster_conversion = NULL,
                    relative_prop = NULL, sample_new_levels = NULL, use_average_cluster_profiles = NULL,
                    rate=1, cluster_size_max=1e24, cluster_size_min=0, sampling_type,
                    nCV = NULL, simulation_parameter = NULL,
                    simulation_model = NULL){
  if (sampling_type=="No_simulation"){
    result=no_simulation(count_table)
    return(result)
  }
  if (sampling_type=="Resampling"){
    result=resampling_from_sc(count_table = count_table, cell_cluster_conversion = cell_cluster_conversion)
    return(result)
  }
  if (sampling_type=="Subsampling"){
    result=subsampling_from_sc(count_table = count_table, rate = rate)
    return(result)
  }
  if (sampling_type=="Subsampling_by_cluster"){
    result=subsampling_by_cluster_from_sc(count_table = count_table, cell_cluster_conversion = cell_cluster_conversion, rate = rate, cluster_size_max = cluster_size_max, cluster_size_min = cluster_size_min, nCV = nCV)
    return(result)
  }
  if (sampling_type=="Simulation"){
    if (simulation_model=="Naive_simulation_with_probe_failure" || simulation_model=="Naive_simulation_without_probe_failure"){
      result=simulation_naive(count_table = count_table, simulation_parameter = simulation_parameter, relative_prop = relative_prop)
    }
    if (simulation_model=="ZINB"){
      result=simulation_ZINB(count_table = count_table, cell_cluster_conversion = cell_cluster_conversion, relative_prop = relative_prop, simulation_parameter = simulation_parameter, sample_new_levels = sample_new_levels, use_average_cluster_profiles = use_average_cluster_profiles)
    }
    return(result)
  }
}


