# Mordred chemical descriptors

A set of ca 1613 chemical descriptors, including both RDKit and original modules. It is comparable to the well known PaDEL-Descriptors (see eos7asg), but has shorter calculation times and can process larger molecules. In this implementation, we fill in empty values using an imputer trained on DrugBank, and eliminate frequently empty columns to a total of 1458 features

This model was incorporated on 2021-09-28.


## Information
### Identifiers
- **Ersilia Identifier:** `eos78ao`
- **Slug:** `mordred`

### Domain
- **Task:** `Representation`
- **Subtask:** `Featurization`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Descriptor`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1458`
- **Output Consistency:** `Fixed`
- **Interpretation:** Vector representation of a molecule

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| abc | float | high | atom-bond connectivity index |
| abcgg | float | high | Graovac-Ghorbani atom-bond connectivity index |
| nacid | float | high | acidic group count |
| nbase | float | high | basic group count |
| spabs_a | float | high | SpAbs of adjacency matrix |
| spmax_a | float | high | SpMax of adjacency matrix |
| spdiam_a | float | high | SpDiam of adjacency matrix |
| spad_a | float | high | SpAD of adjacency matrix |
| spmad_a | float | high | SpMAD of adjacency matrix |
| logee_a | float | high | LogEE of adjacency matrix |

_10 of 1458 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos78ao](https://hub.docker.com/r/ersiliaos/eos78ao)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos78ao.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos78ao.zip)

### Resource Consumption
- **Model Size (Mb):** `319`
- **Environment Size (Mb):** `829`
- **Image Size (Mb):** `1645.45`

**Computational Performance (seconds):**
- 10 inputs: `27.96`
- 100 inputs: `34.52`
- 10000 inputs: `963.71`

### References
- **Source Code**: [https://github.com/mordred-descriptor/mordred](https://github.com/mordred-descriptor/mordred)
- **Publication**: [https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0258-y](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-018-0258-y)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2018`
- **Ersilia Contributor:** [miquelduranfrigola](https://github.com/miquelduranfrigola)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [BSD-3-Clause](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos78ao
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos78ao
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
