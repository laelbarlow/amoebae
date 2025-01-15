
# How to install and run AMOEBAE on Digital Research Alliance of Canada (DRAC) clusters

**Note**: This is mostly relevant for AMOEBAE users in Canada. However, this approach can be easily adapted to any high performance computing cluster.

AMOEBAE depends on Conda packages, but installing Conda packages directly on DRAC isn't appropriate for the system. One way to circumvent that issue is to install dependencies in an Apptainer container, by performing the following steps.

1. Log on to a DRAC cluster login node, clone the amoebae repository, and navigate into the amoebae directory. Upload your sequence data, as needed using an appropriate method (see [DRAC documentation on data transfers](https://docs.alliancecan.ca/wiki/Transferring_data)). Also, I recommend doing the subsequent steps in a screen or tmux session so that your work won't be interrupted if you get disconnected (see [DRAC documentation on prolonging terminal sessions](https://docs.alliancecan.ca/wiki/Prolonging_terminal_sessions)).

2. Build the container using the `pixi.def` definition file provided in the AMOEBAE repository (see [DRAC Apptainer documentation](https://docs.alliancecan.ca/wiki/Apptainer)):
    ```bash
    module load apptainer
    apptainer build --disable-cache pixi.sif pixi.def
    ```
    - This will take a few minutes.

3. To run AMOEBAE workflow steps requiring internet access (`download_queries`
   and `download_dbs`), enter a shell session within the
   container on the login node (unless using the Cedar cluster specifically):
    - Enter the shell session:
        ```bash
        APPTAINER_BIND=''
        apptainer shell -C -B $PWD:/root --pwd /root pixi.sif
        ```
    - Now you are in an environment where all the AMOEBAE dependencies are
      installed. You can run `snakemake` commands as described in the [workflow
      protocol](./workflow_protocol.md), without the need to prefix with `pixi run `.
    - Run the specific rules in the workflow (if/when needed as part of the workflow protocol):
        ```bash
        snakemake --cores 1 download_queries
        snakemake --cores 1 download_dbs
        ```
    - Exit the shell session:
        ```bash
        exit
        ```

4. Otherwise, start an interactive session, so that your analysis will be run on a compute
   node. For example, like this (see [DRAC documentation](https://docs.alliancecan.ca/wiki/Running_jobs#Interactive_jobs)):
    ```bash
    salloc --time=1:0:0 --mem-per-cpu=16G --ntasks=1 --account=def-leppard
    ```

5. Then, enter a shell session within the container:
    ```bash
    module load apptainer
    APPTAINER_BIND=''
    apptainer shell -C -B $PWD:/root --pwd /root pixi.sif
    ```

6. Now run the remaining workflow steps, as described in the workflow protocol. 
    ```bash
    snakemake get_ref_seqs 
    ...
    ```
