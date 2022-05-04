from math import sin, cos
import json
import matplotlib.pyplot as plt
from animation import AxData, Animator


def main():
  with open('swingup_trajectory.json') as f:
    data = json.load(f)

  theta1 = data['q1']
  theta2 = data['q2']

  l1 = 1
  l2 = 1

  pendulum_x = []
  pendulum_y = []

  for i in range(len(theta1)):
    elbowx =  l1 * sin(theta1[i])
    elbowy = -l1 * cos(theta1[i])
    tipx = elbowx + l2 * sin(theta1[i] + theta2[i])
    tipy = elbowy - l2 * cos(theta1[i] + theta2[i])
    pendulum_x.append([0, elbowx, tipx])
    pendulum_y.append([0, elbowy, tipy])


  fig, ax = plt.subplots()
  ax.set_aspect('equal')

  ax_data = [
    AxData(pendulum_x, pendulum_y, 'position', plot_history=False),
  ]

  dt = 0.05
  animator = Animator(1, ax_data, dt, fig, ax, show_legend=False)
  animator.set_save_path('/home/chub/swingup')
  animator.run()

if __name__ == "__main__":
  plt.rcParams['figure.figsize'] = [16, 10]
  plt.rcParams['savefig.facecolor'] = 'black'
  plt.rcParams['figure.facecolor'] = 'black'
  plt.rcParams['figure.edgecolor'] = 'white'
  plt.rcParams['axes.facecolor'] = 'black'
  plt.rcParams['axes.edgecolor'] = 'white'
  plt.rcParams['axes.labelcolor'] = 'white'
  plt.rcParams['axes.titlecolor'] = 'white'
  plt.rcParams['xtick.color'] = 'white'
  plt.rcParams['ytick.color'] = 'white'
  plt.rcParams['text.color'] = 'white'
  plt.rcParams["figure.autolayout"] = True
  # plt.rcParams['legend.facecolor'] = 'white'

  main()

