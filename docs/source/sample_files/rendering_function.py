def orca_render(mendeleev:dict, clusters_cfg:dict, config:dict, file_data:dict, job_specs:dict, misc:dict):
    """Renders the job instructions file and the input file associated with the ORCA program.
    
    Parameters
    ----------
    mendeleev : dict
        Content of AlexGustafsson's Mendeleev Table YAML file (found at https://github.com/AlexGustafsson/molecular-data).
        Unused in this function.
    clusters_cfg : dict
        Content of the YAML clusters configuration file.
    config : dict
        Content of the YAML configuration file.
    file_data : dict
        Information extracted by the scanning function from the geometry file.
    job_specs : dict
        Contains all information related to the job.
    misc : dict
        Contains all the additional variables that did not pertain to the other arguments.
    
    Returns
    -------
    rendered_content : dict
        Dictionary containing the text of all the rendered files in the form of <filename>: <rendered_content>.
    rendered_instructions : str
        Name of the rendered job instructions file, necessary to launch the job.
    
    Notes
    -----
    Pay a particular attention to the render_vars dictionary, it contains all the definitions of the variables appearing in your jinja template.
    """
    
    # Define the names of the templates
    
    template_input = "sample_orca.inp.jinja"
    template_instructions = "sample_orca_job.sh.jinja"
    
    # Define the names of the rendered files
    
    rendered_input = misc['mol_name'] + ".inp"
    rendered_instructions = "orca_job.sh"
    
    # Initialize the dictionary that will be returned by the function
    
    rendered_content = {}
    
    # Render the template for the input file

    print("{:<80}".format("\nRendering the jinja template for the orca input file ..."), end="")
    
    render_vars = {
      "method" : config['method'],
      "basis_set" : config['basis-set'],
      "job_type" : config['job-type'],
      "charge" : config['charge'],
      "multiplicity" : config['multiplicity'],
      "coordinates" : file_data['atomic_coordinates']
    }
    
    rendered_content[rendered_input] = jinja_render(misc['templates_dir'], template_input, render_vars)
    
    print('%12s' % "[ DONE ]")

    # Render the template for the job instructions file

    print("{:<80}".format("\nRendering the jinja template for the orca job instructions file ..."), end="")
    
    render_vars = {
      "mol_name" : misc['mol_name'],
      "config_name" : misc['config_name'],
      "user_email" : config['user-email'],
      "mail_type" : config['mail-type'],
      "job_walltime" : job_specs['walltime'],
      "job_cores" : job_specs['cores'],
      "job_mem_per_cpu" : job_specs['mem_per_cpu'],
      "partition" : job_specs['partition'],
      "set_env" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['set_env'],
      "command" : clusters_cfg[job_specs['cluster_name']]['progs'][job_specs['prog']]['command'],
      "prog" : job_specs['prog']
    }
    
    rendered_content[rendered_instructions] = jinja_render(misc['templates_dir'], template_instructions, render_vars)
    
    print('%12s' % "[ DONE ]")

    # Return the content of the rendered files and the name of the rendered job instructions file
    
    return rendered_content, rendered_instructions