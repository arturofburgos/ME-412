import pandas as pd
import matplotlib.pyplot as plt

# Function that clean the non used CSV part


def df_clean(csv_path):
    df = pd.read_csv(csv_path, header=None).to_numpy()
    max_idx = df[0, :].argmax()
    label = csv_path.split('/')[-1].rstrip('.csv')
    return df[:, 1:max_idx], label


df1, label1 = df_clean(
    'C:/Users/burgo/Desktop/ME-412/#3/data/20 x 20 Jacobi 4.csv')
df2, label2 = df_clean(
    'C:/Users/burgo/Desktop/ME-412/#3/data/20 x 20 Gauss-Seidel 4.csv')
df3, label3 = df_clean(
    'C:/Users/burgo/Desktop/ME-412/#3/data/20 x 20 SOR 4 1_6.csv')
df4, label4 = df_clean(
    'C:/Users/burgo/Desktop/ME-412/#3/data/20 x 20 SOR 4 1_4.csv')
df5, label5 = df_clean(
    'C:/Users/burgo/Desktop/ME-412/#3/data/20 x 20 SOR 4 1_2.csv')

#====================================#
#                                    #
#     Plots and Post-processing      #
#                                    #
#====================================#

fig, ax = plt.subplots()
ax.plot(df1[0, :], df1[1, :], label=label1)
ax.plot(df2[0, :], df2[1, :], label=label2)
ax.plot(df3[0, :], df3[1, :], label=label3)
ax.plot(df4[0, :], df4[1, :], label=label4)
ax.plot(df5[0, :], df5[1, :], label=label5)
ax.set_yscale('log')
ax.set_ylabel('Residual')
ax.set_xlabel('Iteration')
ax.set_title(
    'Comparison between different iterative methods: 20 x 20 and k = 4')
ax.legend()
plt.grid()
plt.show()
