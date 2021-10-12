#!/bin/bash

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/slipstream/home/joeboyd/anaconda2/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/slipstream/home/joeboyd/anaconda2/etc/profile.d/conda.sh" ]; then
        . "/slipstream/home/joeboyd/anaconda2/etc/profile.d/conda.sh"
    else
        export PATH="/slipstream/home/joeboyd/anaconda2/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
exit 0
PATH=/slipstream/home/joeboyd/anaconda2/condabin:$PATH
export PATH
eval $(conda shell.bash hook)

conda init bash
conda activate suppa2_env
