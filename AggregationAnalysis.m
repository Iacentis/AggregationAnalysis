function [AggregationScores,GroupingData,RootMeanSquares] = AggregationAnalysis(Uint8ImageData)
% Uint8ImageData: Image Data, assumed to be in Uint8 format (may work for
% other formats as well, but not tested
% AggregationSquareSize: Size of the square to select aggregation data from
% in pixels, integer
AggregationScores = [0;0];
RootMeanSquares = [0;0];
GroupingData = cell(3,1);
DoubleImageData = double(Uint8ImageData);
Delta = abs(DoubleImageData(3:end,:)-DoubleImageData(1:end-2,:)); %dI / dx
DeltaY = abs(DoubleImageData(:,3:end)-DoubleImageData(:,1:end-2));
Delta = Delta(:,2:end-1)+DeltaY(2:end-1,:);
Delta = Delta./max(Delta(:));
Delta(Delta>mean(Delta(:))) = 1;
Delta(Delta < 1) = 0;
%%%%%%%% Score Calculation

xMeanDeviation = mean(Delta,1);
yMeanDeviation = mean(Delta,2);

GroupingData{1} = Delta;
GroupingData{2} = xMeanDeviation;
GroupingData{3} = yMeanDeviation;

RootMeanSquares(1) = sqrt(sum(xMeanDeviation.^2)./numel(xMeanDeviation));
RootMeanSquares(2) = sqrt(sum(yMeanDeviation.^2)./numel(yMeanDeviation));
AggregationScores(1) = 1./sum(diff([0,find(xMeanDeviation < RootMeanSquares(1))])>1);
AggregationScores(2) = 1./sum(diff([0;find(yMeanDeviation < RootMeanSquares(2))])>1);

%%%%%%% Plot Generation
figure(1)
clf
subplot(2,2,1)
hold off
imagesc(DoubleImageData);
hold on
plot3(1:numel(xMeanDeviation),zeros(size(xMeanDeviation)),xMeanDeviation-min(xMeanDeviation),'LineWidth',2);
plot3(zeros(size(yMeanDeviation))+numel(xMeanDeviation),1:numel(yMeanDeviation),yMeanDeviation-min(yMeanDeviation),'LineWidth',2);
xlabel('X')
ylabel('Y')
axis([0,numel(xMeanDeviation),0,numel(yMeanDeviation),0,max(max(yMeanDeviation)-min(yMeanDeviation),max(xMeanDeviation)-min(xMeanDeviation))])
set(gca,'ztick',[])
view(3)

subplot(2,2,2)
imagesc(smooth(xMeanDeviation,200)'.*smooth(yMeanDeviation,200))

subplot(2,2,3)
hold off
plot(xMeanDeviation);
hold on
plot(xMeanDeviation.*0+RootMeanSquares(1));
xlabel('X')
ylabel('Relative Intensity')
title(sprintf('X-Aggregation Score = %g',AggregationScores(1)))
legend('Data',sprintf('RMS = %g',RootMeanSquares(1)))
axis([0,numel(xMeanDeviation),min(xMeanDeviation),max(xMeanDeviation)])

subplot(2,2,4)
hold off
plot(yMeanDeviation);
hold on
plot(yMeanDeviation.*0+RootMeanSquares(2));
xlabel('Y')
ylabel('Relative Intensity')
title(sprintf('Y-Aggregation Score = %g',AggregationScores(2)))
legend('Data',sprintf('RMS = %g',RootMeanSquares(2)))
axis([0,numel(yMeanDeviation),min(yMeanDeviation),max(yMeanDeviation)])

end