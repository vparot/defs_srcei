function pr = reducepeaks(p)
k1 = [-1 1 0]';
pr = p;
c = 0;
c = c + -conv2(pr,k1,'same'); % positive when next point is higher
[~, ~, th] = isoutlier(diff(p));
c = c - th;
c = ceil(max(0,c/2));
c(end,:) = 0;
pr(~~c) = pr(~~c) + c(~~c);
pr(find(~~c)+1) = pr(find(~~c)+1) - c(~~c);
assert(isequal(sum(p),sum(pr)),'Sum consistency error')