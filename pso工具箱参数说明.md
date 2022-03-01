# pso工具箱

PSO工具箱由美国北卡罗来纳州立大学航天航空与机械系教授Brian Birge开发，名字为[PSOt](D:\jetaime\document\OneDrive - business-cn\学习资料\matlab\MATLAB智能算法30个案例分析\MATLAB智能算法30个案例分析\chapter17 基于PSO工具箱的函数寻优算法\PSOt)

# 工具箱主要函数

| 函数名称               | 函数功能           |
| ---------------------- | ------------------ |
| goplotpso              | 绘图函数           |
| Normmat                | 格式化矩阵数据函数 |
| pso_Trelea_vectorized  | 粒子群优化主函数   |
| linear_dyn，spiral_dyn | 时间计算函数       |
| forcerow.m，forcecol   | 向量转化函数       |

# pso_Trelea_vectorized

`[optOUT，tr，te]=pso_Trelea vectorized（functname，D，mv，VarRange，minmax，PSOparams，plotfcn，PSOseedValue）`

**输入参数：**

1. functname：优化函数名称。
2. D：待优化函数的维数。
3. mv：最大速度取值范围。默认值为4
4. VarRange：粒子群位置取值范围。默认为 -100 到100
5. minmax：寻优参数，决定寻找的是最大化模型、最小化模型还是和某个值最接近。
   当minmax=1时，表示算法寻找最大化目标值；
   当minmax=0时，表示算法寻找最小化目标值；
   当minmax=2时，表示算法寻找的目标值与PSOparams数组中的第12个参数最相近。
6. plotfcn：绘制图像函数，默认为`goplotpso`。
7. PSOseedValue：初始化粒子位置，当PSOparams数组中的第13个参数为0时，该参数有效。
8. PSOparams：算法中具体用到的参数，为一个13维的数组。

**PSOparams中13个参数含义**：

参数默认值为：`PSOparams=[100 2000 24 2 2 0.9 0.4 1500 1e-25 250 NaN 0 0]`

1. PSOparams中的第1个参数表示MATLAB命令窗显示的计算过程的间隔数，100表示算法每迭代100次显示一次运算结果，如取值为0，不显示计算中间过程。
2. PSOparams中的第2个参数表示算法的最大迭代次数，在满足最大迭代次数后算法停止，此处表示最大迭代次数为2000。
3. PSOparams中的第3个参数表示种群中个体数目，种群个体越多，越容易收敛到全局最优解，但算法收敛速度越慢，此处表示种群个体数为24。
4. PSOparams中的第4个参数、第5个参数为算法的加速度参数，分别影响局部最优值和全局最优值，一般采用默认值2。
5. 250
6. PSOparams中的第6个参数、第7个参数表示算法开始和结束时的权值，其他时刻的权值通过线性计算求得，此处表示算法开始时的权值为0.9，算法结束时的权值为0.4。
7.
8. PSOparams中的第8个参数表示当迭代次数超过该值时，权值取PSOparams中的第6个参数和PSOparams中的第7个参数的小值。
9. PSOparams中的第9个参数表示算法终止阈值，当连续两次迭代中对应种群最优值变化小于此阀值时，算法终止，此处值为1e-25。
10. PSOparams中的第10个参数表示用于终止算法的阈值。当连续250次迭代中函数的梯度值仍然没有变化，则退出迭代。
11. PSOparams中的第11个参数表示优化问题是否有约束条件，取NaN时表示为非约束下的优化问题。
12. PSOparams中的第12个参数表示使用粒子群算法类型。
    0 ：普通pso(default)
    1,2 ：Trelea types 1,2
    3 ： Clerc's Constricted PSO, Type 1"
13. PSOparams中的第13个参数表示种群初始化是否采用指定的随机种子，0表示随机产生种子，1表示用户自行产生种子。
