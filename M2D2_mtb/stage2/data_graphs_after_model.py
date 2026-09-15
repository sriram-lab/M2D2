# GRAPHS OF RESULTS
# Bar graph of results for different algorithms
# algorithms = ['Random Forest', 'XGBoost', 'Grad. Boost.']
# r_square = [0.4817, 0.3974, 0.4364]
# std_dev = [0.0964, 0.1035, 0.0844]
# plt.figure(figsize=(8, 5))
# plt.bar(algorithms, r_square, yerr=std_dev)
# plt.title('Pearsons R for Different Algorithms')
# plt.ylabel('Accuracy')
# plt.xlabel('Scores')
# plt.show()

print(f"{np.mean(r_scores):.4f} ± {np.std(r_scores):.4f}")
# LSRL plot of best algorithm predicted vs actual scores 
x, y = y_test, y_pred
plt.scatter(x, y)

# LSRL
A = np.vstack([x, np.ones(len(x))]).T
m, b = np.linalg.lstsq(A, y, rcond=None)[0]
plt.plot(x, m*x + b, 'r',)

plt.xlabel('Actual')
plt.ylabel('Predicted')
plt.title('XGBoost Predicted vs. Actual')
plt.show()