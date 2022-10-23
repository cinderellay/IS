# IS
# <center>System identification experiment



# <center>Single-volume tank system modeling

### 1. The title：single-volume tank system modeling 

### 2. The contents:

The experimental data are collected from the input and output data of the single-volume water tank system, a certain identification algorithm is used to identify the system model, and the results are analyzed.

### 3. The purposes

(1)Learn the theoretical knowledge of system identification and master the method of least squares identification

(2)Analyze and compare different least squares identification methods and their effects

(3)Learn data acquisition and pre-processing during modeling of unknown real systems

(4)Learn to write M files in MATLAB and write least squares identification programs independently

(5)Analyze the deviation between the identified system and the actual system, and verify the experimental results

(6)Add perturbation to the original system , analyze and compare the immunity of different algorithms

### 4. The procedure 

(1)Acquisition of input and output data for single volume water tank systems；

(2)The input and output data obtained are used to identify the order of the known system and estimate its parameters by means of a programmed identification algorithm.The time delay of the system is directly obtained by observing the experimental data, and the correctness of the adopted algorithm is verified by plotting the deviation curve and error variation curve of the identification result and the actual output；

(3)Add perturbation to the system, use the original system parameters to fit the data after adding the perturbation, and observe the anti-perturbation effect

Our group first understands the method of system identification from the theoretical point of view, mastering the principle and algorithm of least squares, and then introducing the recursive least squares method, which greatly reduces the complexity of calculation.Given that the system was not stable at the beginning, it was desired to reduce the importance of the old data and increase the weight of the new data, so the least square method of fading memory was used in this experiment. The method received good experimental results. The fading memory least square algorithm is written, the conclusions are analyzed, and the results are verified. Also write the extended least square method, compare it with the experimental results of the fading memory least squares algorithm, analyze the conclusions, and verify the results. The input and output data that have been obtained are used to estimate the parameters of the unknown water tank system by the identification algorithm. In order to verify the correctness of the experimental results, the deviation curves of the output of the discriminated system and the actual system need to be plotted and the error variation curves need to be plotted to facilitate the observation of the error variation. Finally, perturbation is added to the original system and the perturbation curve is fitted with the parameters of the original system to observe the immunity.


### 5. The principles

#### 5.1 Data preprocessing

The first line of all data files is the time, when importing, the input starts from the second line; to observe the delay from the data, the output data is imported from the subtraction of the delay

#### 5.2 Introduction of least square method

##### 5.2.1 Proposal of linear system model

Considering the system steady-state model,hen the system is running steadily,the relationship between the system input and output at the same moment：

n inputs (independent variables)：$u_i=1,2,……，n$

1 output (dependent variable)：$y$

$N$：Total number of data sets (length of data, number of samples)

For each observation, a set of data is obtained, which constitutes an observation equation. N observations, N sets of data are obtained, which constitutes N observation equations：
$$
y(1)=a_0+a_1u_1(1)+a_2u_2(1)+……+a_nu_n(1)  \\

y(2)=a_0+a_1u_1(2)+a_2u_2(2)+……+a_nu_n(2) \\

  ……  \\

y(N)=a_0+a_1u_1(N)+a_2u_2(N)+……+a_nu_n(N)
$$


Introducing measurement error $ e（k）$​，we get:
$$
y(k)=a_0+a_1u_1(k)+a_2u_2(k)+a_nu_n(k)+e(k)
(k=1,2……N)
$$
$ e（k）$ is also known as observation error, model error, model residuals, noise
$$
e(k)=y(k)-a_0-a_1u_1(k)-a_2u_2(k)-……-a_nu_n(N)\\
\theta=(a_0,a_1,……,a_n)^T\\
Y_N=\Phi_N\theta+e_N
$$
parameter estimates obtained from N data sets：
$$
\hat\theta(N)=(\Phi_N^T\Phi_N)^{-1}\Phi_N^TY_N=\hat\theta
$$
Observation data matrix：
$$
\Phi_N=\left[
\begin{matrix}
\phi_1^T\\
\cdots\\
\phi_{n-1}^T \\
\phi_n^T
\end{matrix}
\right]
=\left[
\begin{matrix}
1 & \cdots & u_{n-1}(1)& u_n(1)\\
\cdots & \cdots & \cdots & \cdots\\
1 & \cdots & u_{n-1}(N-1) &u_n(N-1)\\
1 & \cdots & u_n(N) & u_n(N)
\end{matrix}
\right]
$$

##### 5.2.2 Recursive least square method

**Introduction to recursive least square method**

Least-squares identification of linear steady-state models using online data, with observation equations established from N sets of old data：
$$
Y_N=\Phi_N\theta+e_N
$$
Parameter estimates obtained from N sets of old data：
$$
\hat\theta(N)=(\Phi_N^T\Phi_N)^{-1}\Phi_N^TY_N=\hat\theta
$$
New set of observations obtained：
$$
y(N+1)\\
u_1(N+1),……,u_n(N+1)
$$
We let$P_N=(\Phi^T_N\Phi_N)^{-1}$

Then$P_{N+1}=(P_N^{-1}+\Phi_{N+1}\Phi_{N+1}^T)^{-1}$

We finally obtain $n+1 $  parameter estimates， 
$$
\theta(N+1)=\theta(N)+\frac{P_N\phi_{N+1}\phi_{N+1}^T}{1+\phi_{N+1}^TP_N\phi_{N+1}}[y(N+1)-\phi_{N+1}^T\phi_{N+1}^{T}\theta(N)]
$$
we call $[y(N+1)-\phi_{N+1}^T\phi_{N+1}^T\theta(N)]$ **innovation**，which is the fitting error to the new observations when using the old estimate $\theta(N)$,the new estimate of the parameter $\theta$ = the old estimate + the corrected value。

We let
$$
K_{N+1}=\frac{P_N\phi_{N+1}\phi_{N+1}^T}{1+\phi_{N+1}^TP_N\phi_{N+1}}\\
\theta(N+1)=\theta(N)+K_{N+1}[y(N+1)-\phi_{N+1}^{T}P_N]\\
P_{N+1}=P_N-K_{N+1}\phi_{N+1}^TP_N
$$
**Selection of estimated initial values $\hat\theta(0)$ and $P_0$**

1、Apply batch least squares to the first N sets of data:
$$
\hat\theta(N)=(\Phi_N^T\Phi_N)^{-1}\Phi_N^TY_N
$$

$$
P_N=(\Phi_N^T\Phi_N)
$$

 Use this as the initial value $\hat\theta(0)$，$P_0$，and recursive from step N+1

2、Select：$\hat\theta(0)=(0,0,……，0)^T$,$P_0=\alpha I$,$\alpha>0$，and consider $\alpha$ to be a sufficiently large number, generally taken as $10^2-10^6$

Taking this as the initial value, after $N$ steps of recursion, we have:
$$
P_N^*=[P_0^{-1}+\Phi_N^T\Phi_N]^{-1}
$$
Notice that$\lim\limits_{\alpha\rightarrow\infty}P_0^{-1}=\lim\limits_{\alpha\rightarrow\infty}\frac{1}{\alpha}I=0$

i.e.$\lim\limits_{\alpha\rightarrow\infty}P_N^{*}=[\Phi_N^T\Phi_N]^{-1}=P_N$

When $\alpha$ is large enough, the initial value $P_0=\alpha I$ has little effect on $P_N$
$$
\hat\theta(N)=P_N^*[\Phi_N^TY_N+P_0^{-1}\hat\theta_0]
$$

$$
\because\lim\limits_{\alpha\rightarrow\infty}P_0^{-1}=0
$$

$$
\lim\limits_{\alpha\rightarrow\infty}\theta(N)=\lim\limits_{\alpha\rightarrow\infty}P^*[\Phi_N^T\Phi_N]=P_N\Phi_N^TY_N=(\Phi_N^T\Phi_N)^{-1}\Phi_N^TY_N
$$

When $\alpha$ is large enough, the initial values $P_0=\alpha I$, $\hat\theta_0=0$ have little effect on $\hat\theta(N)$, there is no need to worry about the specific values, the algorithm itself is convergent

##### 5.2.3 Fading memory least squares

As the data grows, $P(k)$ becomes smaller and smaller and finally tends to zero, and the least squares recognition algorithm gradually loses its correction capability, which is a phenomenon called data saturation. Because the $P(k)$ matrix is used to accumulate the amount of information, as time grows, the information of old data accumulates in the matrix $P(k)$, and at a certain point the information of new data cannot be added in, so the algorithm will gradually lose its correction ability. For time-varying systems, due to the accumulation of old data information, the least-squares recognition algorithm will also gradually lose its tracking ability for the same reason. Both the data saturation problem and the time-varying tracking capability problem are due to the accumulation of data information, which is not decayed and the amount of information obtained from new data decreases, affecting the correction capability and tracking capability of the algorithm.  The forgetting factor method is an identification method proposed to overcome the data saturation phenomenon and solve the time-varying tracking problem,and its basic idea is to add forgetting factors to the old data to reduce the proportion of old data information in the matrix $P(k)$ and increase the content of new data information. For $N$ sets of observations, there are objective functions：
$$
min\quad J_N(\theta)=\sum_{k=1}^{N}\lambda ^{N-k}e^2(k)
=e^2(N)+\lambda e^2(N-1)+……+\lambda ^{N-1}e^2(1)\\
$$
where $\lambda$ is the forgetting factor, $0<\lambda \leqslant 1$,In the new objective function, due to the introduction of the forgetting factor less than 1, the fitting error to the new data is minimized to be emphasized, while the fitting error to the old data is gradually faded away。

Then we get the final fading memory least squares identification algorithm ：
$$
\hat\theta _{N+1}=\hat\theta_N+K_{N+1}[y(N+1)-\phi_{N+1}^{T}\hat\theta_N]\\
K_{N+1}=\frac{P_N\phi_{N+1}}{\lambda +\phi_{N+1}^{T}P_N\phi_{N+1}}\\
P_{N+1}=\frac{1}{\lambda}[P_N-K_{N+1}\phi_{N+1}^TP_N]
$$
Forgetting factor $0<\lambda \leqslant 1$, generally take $0.95\leqslant \lambda \leqslant 1$, when $\lambda$=1,it becomes recursive least square. The smaller the forgetting factor, the more importance is attached to new data, the heavier the weighting of new fitting errors, the better the algorithm can track the changes of system parameters; the smaller the forgetting factor, the greater the influence of parameter estimation receiving observation noise, interference, and thus the greater the fluctuation of parameter estimates

##### 5.2.4 Extended least squares

Extended least square method：The least square parameter estimation algorithm is unbiased when the noise mean is 0. In order to solve the biased nature of least square parameter estimation, the identification of the noise model is taken into account at the same time, which is the extended least square method. The extended least square method is equivalent to the least square method with expanded dimensionality of the parameter vector and data vector, which can give consistent estimates of the parameters in the case of colored noise (represented by the mean sliding model), and at the same time, the noise model can also be identified。

If the error is measurable, the error can be expanded into the input vector, and the expanded vector is called the augmented vector:
$$
\phi_{t-1}=[-y_{t-1},……,-y_{t-p},\epsilon_{t-1},……,\epsilon_{t-q}]^T\\
\theta=[a_1,……,a_p,b_1,……,b_q]^T
$$
Then we have $y_t=\phi_{t-1}^T\theta+\epsilon_t$ 。The recursive least square method can be used to solve for the coefficients θ。But in practice, the error $ \epsilon$ at each step is not measured exactly and $ \phi_t$ is not exact, and thus can only be replaced by the estimated amount of $ \phi$：
$$
\hat\phi_{t-1}=[-y_{t-1},……,-y_{t-p},\hat\epsilon_{t-1},……,\hat\epsilon_{t-q}]^T
$$
Thus the estimator $\hat{\theta}_t $ is derived.Usurally $\hat{\theta}_t $ converges to $\theta$, so the method is valuable！

The fitting error $e_t$ at moment t is used as an estimate of the unknown unpredictable noise $\xi_t$：$e_t=y_t-\phi_t^T\hat\theta_{t-1}$（预测形式）

Parameter Variable：$\theta=(a_1,a_2,……,a_{n_a},b_0,b_1,……,b_{n_b},c_1,c_2,……，c_{n_c})^T$​

The estimated parameter vector at moment t:$\hat\theta=(\hat a_1,\hat a_2,……,\hat a_{n_a},\hat b_0,\hat b_1,……,\hat b_{n_b},\hat c_1,\hat c_2,……，\hat c_{n_c})^T$​

The augmented observation vector at time t：$\phi_t=(-y_{t-1},……,-y_{t-n_a},u_{t-k},……,u_{t-k-n_b},e_{t-1},……,e_{t-n_c})^T$

Estimated value of unmeasurable noise at moment t$e_t$:
$$
e_t=
\begin{cases}
0,\quad x< 0\\
y_t-\phi_t^T\hat \theta_t, \quad x\geq 0
\end{cases}\\
K_t=\frac{P_{t-1}\phi_t}{\lambda +\phi_t^TP_{t-1}\phi_t}\\
\hat\theta_t=\hat\theta_{t-1}+K_t(y_t-\phi_T^T\hat \theta_{t-1})\\
P_t=\frac{1}{\lambda}(I-K_t\phi_t^T)P_{t-1}
$$

#### 5.3 Design of the algorithm

##### 5.3.1 Fading memory least squares

In designing the fading memory least square method, we use the $CAR$ model and choose a second-order system
$$
CAR 　Mode:y(k)=-\sum_{i=1}^{n}a_iy(k-1)+\sum_{i=1}^{n}b_iu(k-i)+w(k)\\
Y_N=\Phi\theta+e_N
$$
For second-order systems, the recursive least squares of fading memory：parameter vector $\theta=[a_1,a_2,b_1,b_2]$,observation vector：$\Phi(k)=[-y(k-1),-y(k-2),u(k-1),u(k-2)]^T$

In this experiment, we add the forgetting factor $\lambda$,and get the new matrix：
$$
\theta(k)=\theta(k-1)+K(k)[y(k)-\phi^T(k)\theta(k-1)]\\
K(k)=P(k-1)\phi(k)[\phi^T(k)P(k-1)\phi(k)+\lambda]^{-1}\\
P(k)=[I-K(k)\phi^T]P(k-1)/\lambda
$$
If $\lambda =1$, it becomes recursive least square。 

##### 5.3.2 Extended least squares

Design the incremental least squares algorithm and choose the controlled autoregressive sliding average $CARMA$ model for the second order system
$$
CARMA　 Mode:y(k)=-\sum_{i=1}^{n}a_iy(k-1)+\sum_{i=1}^{n}b_iu(k-i)+\sum_{i=1}^{n}c_i\xi(k-i)+w(k)\\
$$
For second-order systems，Extended least square：parameter vector： $\theta=[a_1,a_2,b_1,b_2,c_1,c_2]$

observation vector：$\Phi(k)=[-y(k-1),-y(k-2),u(k-1),u(k-2),\xi(k-1),\xi(k-2)]^T$

$RELS$ identification algorithm:
$$
\hat\theta(k)=\hat\theta(k-1)+K(k)[z(k)-h^T(k)\hat\theta(k-1)]\\
K(k)=P(k-1)h(k)[h^T(k)P(k-1)h(k)+1]^{-1}\\
P(k)=[I-K(k)h^T]P(k-1)
$$
The initial value is taken as：$P_0=aI$，$a$ is taken as $10^2-10^6$，$\hat\theta_0=\theta_0=(0　0 　\cdots　0)^T$

### 6. The results and analysis

####  6.1 Identification results

##### 6.1.1 Fading memory recursive least squares 

We choose the average of the last five sets of data of the identification parameters as the system model parameters

|      参数       |  $a_1$  |  $a_2$  | $b_1$  | $b_2$  |
| :-------------: | :-----: | :-----: | :----: | :----: |
|   $\lambda=1$   | -0.8101 | -0.1458 | 0.0032 | 0.0032 |
| $\lambda=0.996$ | -0.7865 | -0.1688 | 0.0035 | 0.0035 |

<figure class="half">
    <img src="C:\Users\AB158\Desktop\图片\递推参数.png"width=290>
    <img src="C:\Users\AB158\Desktop\图片\渐消参数.png"width=290>
</figure>


In the design of the program, the selection of the matrix P must be positive definite and larger values, and after several tests, it was found that 10^6 was the best result, which is also consistent with the theoretical experience. In the initialization of $\theta$ parameters, we first chose the initial value of all 3, and found that the effect is not very good, and the number of iterations is not enough when there is not enough amount of data, which will lead to serious deviation of the parameters, and finally we chose 0 as the initial value. Because in this experiment, the delay is obtained by observation of the data, but there is a certain error in the delay of each group of experimental data, which may be caused by human control at the beginning of the experimental stage, so in the debugging process, the best effect of this experiment is around 20s, of course, there are data with a relatively small delay, for example, in the subsequent verification, there is a group of test data with a delay of about 8s. The wrong selection of time delay will lead to the inclusion of unnecessary interference data during training, which will affect the recursive effect of the data, and will lead to serious deviation of the discriminative data from the real data in the data fitting stage, so the time delay has a great impact on the effect of the experiment. Finally, comparing the recursive least squares and the fading memory least squares, the difference between the identified parameters is not large, and we choose the fading memory algorithm that can play a greater role in the new data as the subsequent validation parameters. The error is $e^2=0.0025878040$ Variance :$0.0000066868$

##### 6.1.2 Extended least squares

We choose the average of the last five sets of data of the identification parameters as the system model parameters

|  $a_1$  |  $a_2$  | $b_2$  | $b_2$  |  $c_1$  |  $c_2$  |
| :-----: | :-----: | :----: | :----: | :-----: | :-----: |
| -0.7911 | -0.1656 | 0.0034 | 0.0034 | -0.0052 | -0.0060 |

<figure class="half">
    <img src="C:\Users\AB158\Desktop\图片\增广参数.png"width=580>
</figure>

In designing the extended least squares, it is similar to the recursive least algorithm, except that a white noise model is added to the model of the recursive algorithm, and the iterative formula is modified to include noise in the parameter calculation. In the subsequent validation, a set of data parameters with smaller error in the experiment is selected as the subsequent validation parameters. Augmentation algorithm error: $e^2=0.0065536480$ variance: $0.0000171114$

#### 6.2 Result verification

##### 6.2.1 Fading memory least squares

The selected parameters：$\theta=[-0.7865,-0.1688,0.0035,0.0035]^T$

The system model：$y(k)=-0.7865y(k-1)-0.1688y(k-2)+0.0035u(k-1)+0.0035u(k-2)$

The model was used to identify the original system and the following comparison chart was obtained

<figure class="half">
    <img src="C:\Users\AB158\Desktop\图片\渐消验证.png"width=290>
    <img src="C:\Users\AB158\Desktop\图片\渐消验证误差.png"width=290>
</figure>


​         												error $e^2=0.2136594985$   variance $0.0004460532$

<figure class="half">
    <img src="C:\Users\AB158\Desktop\图片\渐消验证误差2.png"width=290>
    <img src="C:\Users\AB158\Desktop\图片\渐消验证2.png"width=290>
</figure>


​						 								error $e^2=0.0710323991$   variance $0.0002006565$

Comparing the two sets of tests, the data fitted by changing the input and affecting the output to produce step changes during the experiment is not as good as the one without changing the input, but the errors of both sets of tests are not large, the final error of the first set of data fluctuates around 0.02, and the error of the second set of data eventually tends to 0, which can indicate that the identified model parameters are reasonable

##### 6.2.2 Extended Recursive Least Squares

The selected training parameters:$\theta=[-0.7911,-0.1656,0.0034,0.0034,0.0052,0.0060]^T$

The system model：$y(k)=-0.7911y(k-1)-0.1656y(k-2)+0.0034u(k-1)+0.0034u(k-2)+0.0052\xi(k-1)+0.0060\xi(k-2)$

<figure class="half">
    <img src="C:\Users\AB158\Desktop\图片\增广验证2.png"width=290>
    <img src="C:\Users\AB158\Desktop\图片\增广验证误差2.png"width=290>
</figure>


​	 													error $e^2=0.1757348177$   variance $0.0003668785$

<figure class="half">
    <img src="C:\Users\AB158\Desktop\图片\增广验证.png"width=290>
    <img src="C:\Users\AB158\Desktop\图片\增广验证误差.png"width=290>
</figure>


​	 													error $e^2=0.0324479973$   variance $0.0000916610$

The effect of using the parameters identified by the extended recursive algorithm for verification is not much different from the fading memory, and the trend is more or less the same, and also fits the original system model better, and the error is also within a certain tolerance range, but by calculating the error we can find a small difference, the error obtained by the extended recursive algorithm is a little smaller, and the results obtained from multiple experiments are also the same, which can show that the generalization ability of the extended recursive algorithm is relatively better

#### 6.3 Verify the data added to the perturbation

In this experiment, we choose to slightly rotate the valve after the second change of input (i.e., at the second step) and then quickly restore it to obtain the raised perturbation as shown in the figure with red arrows, because the two algorithms used in this experiment, fading memory recursion and extended recursion, both belong to recursive algorithms and the parameters are updated in real time, so they have strong immunity to perturbation, a feature that can also be found in the following verification.

The same parameters were taken as for the normal data test with no interference as described above.

Fading memory model：$y(k)=-0.7865y(k-1)-0.1688y(k-2)+0.0035u(k-1)+0.0035u(k-2)$

Augmented model：$y(k)=-0.7911y(k-1)-0.1656y(k-2)+0.0034u(k-1)+0.0034u(k-2)+0.0052\xi(k-1)+0.0060\xi(k-2)$

<figure class="half">
    <img src="C:\Users\AB158\Desktop\图片\渐消.png"width=290>
    <img src="C:\Users\AB158\Desktop\图片\渐消误差.png"width=290>
</figure>


​														error$e^2=0.2704154898$   variance $0.0006676926$

<figure class="half">
    <img src="C:\Users\AB158\Desktop\图片\增广.png"width=290>
    <img src="C:\Users\AB158\Desktop\图片\增广误差.png"width=290>
</figure>


​														error$e^2=0.2329327028$   variance $0.0005751425$

Compared with the unperturbed data, it can be seen from the error that the result is obviously worse after adding the perturbation, which is reasonable because it cannot be fitted properly at the perturbation, which will certainly lead to an increase in the error, but according to the validation images and the calculated error, it can also be found that there is still a better identification effect for the perturbed data, and the final error also tends to 0. Only at the perturbation there is a larger error change, which can show that the recursive algorithm has strong anti-perturbation ability. Comparing the two algorithms , the fit of the two algorithms is not much different and the error variation is similar, but the fading memory recursive curve is fitted more perfectly, however, the starting error is larger and the final error is also relatively larger. In order to have a more obvious comparison, a second set of data was used again for validation (here the validation graph is not added, and the data is named: rd_2.mat) to obtain the extended recursion: error:$ 0.4226276311$ variance:$ 0.0008860118$; the fading memory recursion: error $ 0.4878451336$ variance $ 0.0010227361 $, the results show that the error of the extended recursion is still smaller, which can indicate that the extended recursive least square algorithm has better immunity to perturbation, i.e., more robustness, compared to the fading memory recursive least square algorithm.

### 7. The summary

Thanks to Mrs. Qin's kind guidance in this experiment, our group has gained a lot. Through this experiment, we have a better understanding of system identification. System identification is a theory and method to study how to use the experimental or operational input and output data of a system containing noise to build a mathematical model of the object under study. System identification is a method of modeling, and different subject areas correspond to different mathematical models. In a sense, the process of development of different disciplines is the process of building the mathematical model. The problem of identification can be reduced to an algorithm of representing the essential characteristics of an objective system (or a system to be constructed) by a model, and using this model to represent the understanding of the objective system into a useful form. There are three elements of identification: data, model, and criterion. Identification is the selection of a model that fits the data better among a set of model classes according to a criterion. In short, the essence of identification is to select a model from a set of model classes that best fits the static or dynamic properties of the actual process of interest according to some criterion. In this experiment, Mr. Qin gave us the specific steps of the experiment in detail before the experiment and kept answering our questions during the process of our experiment. There are some flaws in our experimental process. First, we are very arbitrary in the process of collecting data, simply and crudely thought that the data are arbitrarily more representative, collected a lot of strange data, and did not pay attention to the sampling period, the sampling period is not controlled to the same value, and in the increase of perturbation, the data did not reach stability to arbitrarily add interference, resulting in very low data differentiation. Therefore, at the initial design stage of the algorithm, the training was able to obtain parameters with low errors, but when the validation was carried out, it always failed to fit the original system better and the errors were particularly large, which was not solved for a long time. In the discussion with classmates we also did not find problems with the algorithm, and eventually chose to use their data, only to find that there was a problem with the original data we collected, so we rescheduled an appointment with the teacher and recollected the data, at which point we were able to complete the basic experimental requirements. Second, we were not very familiar with various least square algorithms, which led to the slow progress of the experiment, and we did not have a very thorough understanding of the algorithm principles, and spent a lot of time in debugging to get better results. Third, when we did the experiment for the first time, we did not recognize the basic principle of least square algorithm for system identification, we initially explored in continuous system, and were confused to do the experiment until the second time when we collected data to really understand that the algorithm for identification is discrete difference equation. So before conducting an experiment, we should learn more clearly the basic principles of the algorithm used in the experiment in order to better help complete the experiment. Our group gave encouragement to each other during the learning process and worked together as a team to finally complete the experimental operation and code refinement. Through this experiment, our group gained a deeper understanding of both system identification and model verification. Once again, we sincerely thank Mr. Qin for his teaching!

### 8. The appendix

~~~ matlab
Appendix1
%MLS.m
%Fading memory recursive least squares 
clc,clear;
close all;
load train.mat
datalength=data(1,1);%total data length
delay=20;%time delay
mu=0.996;%forgetting factor
length=datalength-delay;%total length of data after subtracting the time delay
u=data(2:length+1,2);%input
y=data(2+delay:datalength+1,3);%output

%Solving parameters using recursive algorithms
%Parameter initialization
theta = zeros(4,1);
Theta = zeros(4,1);
J=0;
P = 10^6*eye(4);%P is generally chosen as 10^6
for k = 3:length-1
    h = [-y(k-1) -y(k-2) u(k-1) u(k-2)]';%Iterate through each input and output import of the second-order system
    K = P*h/(mu + h'*P*h);      %updatie K
    P = (P - K*h'*P)/mu;        %updatie P
    J=J+(y(k) - h'*theta)^2/(mu + h'*P*h);  %calculate the criterion function
    theta = theta + K*(y(k) - h'*theta);    %update identification parameters:a1,a2,b1,b2
    Theta = [Theta,theta];                  %storage identification parameters
end
%%Calculating errors and plotting parameter curves
L=length-2;%true data length
y1=y(3:length);%true value
y2=zeros(L,1);%estimated value
a=size(y2);
for i=1:2 
 y2=y2-Theta(i,a(1))*y(3-i:length-i);%calculate the output identification tree pool first
end
for i=1:2 
 y2=y2+Theta(i+2,a(1))*u(3-i:length-i);%calculate the discriminative output of the input
end
error=y1-y2;%residual vector
J=error'*error;
J2=J/L;%residual variance
if(mu==1)
     fprintf("递推最小二乘法 RLS\n");
else fprintf("渐消记忆递推最小二乘法 MLS\n");
end
fprintf("误差：J = %.10f\n",J);
fprintf("方差：J2 = %.10f\n",J2);

sizeiden=size(Theta);
k=1:sizeiden(2);
plot(k,Theta');grid;
if(mu==1)
     title("加扰动递推最小二乘参数变化 ");
else title('加扰动渐消记忆递推最小二乘法')
end
legend('a1','a2','b1','b2');

~~~

~~~matlab
附录2
%ELS.m
%Extended least square recursive algorithm 
clear;clc;
close all
%% Importing data
load train.mat
datalength=data(1,1);%total data length
u=data(2:datalength,2);%data input
v=randn(1,datalength);%construct random perturbations
delay=20;%time delay
length=datalength-delay;%total length of data after subtracting the time delay
y=data(2+delay:datalength,3);

%% Recursive solution
P=10^6*eye(6); %the initial value of P is chosen between 100 and 10^6, generally 10^6 is chosen
Theta=zeros(6,length-2);     %estimated values of the parameters, storing the intermediate process valuation
% K=[10;10;10;10;10;10];     %initialize the K matrix
for i=3:length-1
    h=[-y(i-1);-y(i-2);u(i-1);u(i-2);v(i-1);v(i-2)];%iterate through each input, output, and perturbation import of the second-order system
    K=P*h*inv(h'*P*h+1);    %update K
    P=(eye(6)-K*h')*P;      %update P
    Theta(:,i-1)=Theta(:,i-2)+K*(y(i)-h'*Theta(: ,i-2));%calculate parameter vectors[a1,a2,b1,b2,c1,c2]，all in the Theta matrix
    v(i)=y(i)-h'*Theta(:,i-1);  %update Perturbation
end

%% Calculate the error and plot the parameter identification results
L=length-2;%true data length
y1=y(2:length-1);%true value
y2=zeros(L,1);%estimated value
a=size(y2);
for i=1:2 %1-2
 y2=y2-Theta(i,a(1))*y(3-i:length-i);%first calculate the first two orders of the output
end 
for i=1:2 
 y2=y2+Theta(i+2,a(1))*u(3-i:length-i);%calculate the first two orders of the input
end
error=y1-y2;%residual vector
J=error'*error;
J2=J/L;%residual variance
fprintf("增广最小二乘法\n");
fprintf("误差：J=%.10f\n",J);
fprintf("方差：J2=%.10f\n",J2);

%Plot identification parameter curves
sizeiden=size(Theta);
i=1:sizeiden(2);
figure (1)
plot(i, Theta(1, :), i,Theta(2,:) , i,Theta(3,:), i, Theta(4,:), i, Theta(5,:), i,Theta(6,:))
title('增广最小二乘法')
legend('a1','a2','b1','b2','c1','c2');
~~~

~~~matlab
附录3
%test.m
%Effect validation
close all
clearvars -except Theta   %Only the identified parameters are retained for testing

%% The average of the last five data of the parameter is taken as the test parameter
sizeiden=size(Theta);
ls=Theta(1:4,sizeiden(2)-4:sizeiden(2));
ls=mean(ls,2);

%% Read and process data
load rd_2.mat
datalength=data(1,1); %total data length
delay=18;  %observe the obtained delay time, test1 delay is 20, test2 delay is 8

length=datalength-delay-2;  %length of real output data
u=data(2:datalength+1,2);%input
y=data(2:datalength+1,3); %output
u0=u(3:length+2);%u(t)
u1=u(2:length+1);%u(t-1)
u2=u(1:length);%u(t-2)

%% Fit the original data to get the identified output
y2(1:delay+2,1)=y(1:delay+2);
for k=3:datalength
    i=k-delay-2;  %the fit is guaranteed to start from the time after the delay
    if(i>0)       
        y2(k,1)=[-y2(k-1,1) -y2(k-2,1) u1(i) u2(i)]*ls;
    else          %the input before the time delay is set as the initial input
        y2(k,1)=[-y2(k-1,1) -y2(k-2,1) 6.5 6.5]*ls;   
    end
end
%% Comparison of the results of the original identification system
t=1:datalength;
figure(1);
plot(t,y);hold on;
plot(t,y2);hold off;
grid;
title("渐消记忆递推验证对比");
legend('实际输出','辨识系统的输出');

figure(2)
yerror=y-y2;    %calculate the error of the identification curve and the original curve
plot(t,yerror);grid; J=yerror'*yerror;%test error J=e^2
Jvc=J/datalength;   
title("渐消记忆递推辨识系统与实际系统误差曲线 e");    
fprintf("检验误差：J = %.10f\n",J);
fprintf("检验方差：Jvc = %.10f\n",Jvc);
~~~

