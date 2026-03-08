【？】immune_checkpoint分析目前绘图时，如果分组匹配不完全，绘图会出错

# 一、关于前端页面

## 五种分析都需要用户上传1个或者2个csv：
框1：数据表 CSV（必传）。
框2：分组表 CSV（选传，不填就自动识别）。
上传分组文件的地方可以加一行提示：“请上传 CSV 格式，第一列为样本名，第二列为分组信息（如 0/Control/Normal/healthy/1/Tumor/Cancer/Treat等）

接下来开始配置，进入LINUX命令行，按照以下提示输入命令即可

# 二、miniconda，目标：用于管理R环境

## 1.下载安装包
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

## 2.执行安装 (遇到enter就回车，遇到提示[yes/no]输入 yes 并按回车)
bash Miniconda3-latest-Linux-x86_64.sh

## 3.刷新环境变量，或者重启终端也行
source ~/.bashrc

命令行前出现(base)则成功


## 4.配置国内镜像（四行命令）

conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/

conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/

conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/

conda config --set show_channel_urls yes


## 5.建立环境

### 创建环境，这里环境名是r_doda

conda create -n r_doda r-base=4.3.3 bioconductor-edger bioconductor-limma -y

### 激活环境
conda activate r_doda

命令行前的（）会显示你所在环境，例如（r_doda）

### 运行完上面俩指令即可进入三，后面这些仅供查阅无需运行

退出当前环境
conda deactivate

查看所有环境
conda env list

删除环境
conda remove -n (环境名) --all

包管理
        命令	              说明               举例	
conda install 包名	     安装包（当前环境） conda install numpy
conda install 包名=版本	  安装指定版本包	 conda install tensorflow=2.6.0
conda list	         列出当前环境已安装的包	 conda list
conda search 包名	      搜索可用包版本	conda search pandas
conda update 包名	       更新某个包	    conda update numpy
conda update --all	       更新所有包	    conda update --all
conda remove 包名	       卸载某个包	    conda remove numpy



# 三、配置R环境(miniconda)

## 安装其他需要的包（limma和edgeR前面已经附带了）

## 安装前确保命令行前有(r_doda)，这样才会把包装到这个环境里面

conda install bioconductor-edger bioconductor-limma bioconductor-deseq2 -y

conda install r-ggpubr r-car r-lme4 r-nloptr bioconductor-enhancedvolcano -y

conda install bioconductor-clusterprofiler bioconductor-org.hs.eg.db -y

conda install bioconductor-gseabase bioconductor-gsva -y

# 四、目前目录结构

doda/
├── run_deg.sh               # 差异分析入口脚本
├── ....sh                   # 其他分析
├── Rscripts/                # 内部 R 脚本逻辑 (无需调用)
├── uploads/                 # 用户上传区 (按 userId 分区)
│   └── {userId}/            
│       ├── data.csv
│       └── group.csv
└── results/                 # 结果输出区
    └── {userId}/{taskID}/   
        ├── analysis.log     # 运行日志 (分析出错可查阅)
        ├── result.csv       # 数据结果
        └── plot.png         # 导出的图表

# 分析测试命令（test_group是test_real_data的分组表，这俩是一组；test_matrix是另一个测试数据）
# DEG吃前者，cibersort和immune_checkpoint吃后者，GSVA和estimate两者都能分析

1.DEG
bash run_deg.sh 01 test_data/test_real_data.csv

2.estimate
bash run_estimate.sh 01 test_data/test_real_data.csv
bash run_estimate.sh 03 test_data/test_matrix.csv

3.GSVA
bash run_GSVA.sh 01 test_data/test_real_data.csv
bash run_GSVA.sh 03 test_data/test_matrix.csv

4.immune_checkpoint
bash run_ICG.sh 03 test_data/test_matrix.csv

5.cibersort
bash run_cibersort.sh 03 test_data/test_matrix.csv

# 给脚本加执行权限（个人本地测试可忽略）

先进入doda目录

赋予入口BASH脚本执行权限
chmod +x *.sh

赋予数据目录最高读写权限 (防止 Java 写入失败)
chmod -R 777 uploads/
chmod -R 777 results/
