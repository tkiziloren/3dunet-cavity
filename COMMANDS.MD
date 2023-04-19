# Creating Conda Environment.

``` bash
conda env create -n cavity --file env.yaml 
```

# Updating existing environment(after updating/changing env.yaml).
``` bash
conda env update --file env.yaml --prune
```