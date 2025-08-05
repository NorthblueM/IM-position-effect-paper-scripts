
## 使用方法

* 运行
    * `cmd`进入到`imide_position`文件夹
    * 命令：`python flow.py`

* 参数路径
    * 都在`flow.py`的`_main()`函数中

### 氨基酸亚胺离子
* 主程序：`flow.py`
* 需求
    * 一个序列中只含一个该氨基酸
* 参数设置
    * `modi_one = set([])`: `modi_one`为空，即不是统计修饰的亚胺离子模式
        * 修饰的亚胺离子模式，设置`aa_one_ls = []`
    * `aa_one_ls`: 要统计的氨基酸种类列表
    * `dct_chrc`: 各氨基酸亚胺离子质量（不带电荷）
    * `aa2chrc`: 氨基酸和亚胺离子符号对应关系（对照`aa_one_ls`,`dct_chrc`）
    * `pep_level`: 最后结果`_label.spectra`文件是否是肽段层次
        * 肽段层次：一个肽段保留排名最靠前的第一个PSM
        * **建议开启，有的肽段采集了很多PSM，会影响位置效应判定**
* 结果
    * 最终的标注结果文件`_label_X.spectra`, `X`是指代具体某氨基酸


### 多能量触发使用方法
* pFind搜索时需要导出mgf
* 谱图文件需要导出一些头信息
    * 运行: `python mgf_file_df.py`
    * 参数: 修改`_main()`函数中存放`pf2`文件夹的路径`dpath`即可（应该与raw同目录），注意是**文件夹**
* 数据处理
    * 运行: `python flow_multi_energy.py`
    * 参数：都在`flow_multi_energy.py`的`_main()`函数中
* 画图
    * `plot_multi_energy.ipynb`
    * 参数：都在jupyter notebook文件的开头

