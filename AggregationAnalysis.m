function [AggregationScores,Plots,GroupingData] = AggregationAnalysis(Uint8ImageData,AggregationSquareSize)
% Uint8ImageData: Image Data, assumed to be in Uint8 format (may work for
% other formats as well, but not tested
% AggregationSquareSize: Size of the square to select aggregation data from
% in pixels, integer

DoubleImageData = double(Uint8ImageData);
%DoubleImageData = DoubleImageData./max(DoubleImageData(:));
Delta = DoubleImageData(2:end,:)-DoubleImageData(1:end-1,:); %dI / dx
Delta = Delta(:,2:end)-Delta(:,1:end-1); % dI^2 / dx dy
SquareData = zeros(size(Delta,1),size(Delta,2),(AggregationSquareSize*2+1)^2);
I = 1;
for Range1 = -AggregationSquareSize:AggregationSquareSize
    for Range2 = -AggregationSquareSize:AggregationSquareSize
        SquareData(:,:,I) = circshift(Delta,[Range1,Range2]);
        I = I+1;
    end
end
DeviationData = std(SquareData,0,3);

%%%%%%%% Score Calculation

xMeanDeviation = mean(DeviationData,1);
yMeanDeviation = mean(DeviationData,2);

xMeanDeviation = xMeanDeviation-min(xMeanDeviation(:));
xMeanDeviation = xMeanDeviation./max(xMeanDeviation);
yMeanDeviation = yMeanDeviation-min(yMeanDeviation);
yMeanDeviation = yMeanDeviation./max(yMeanDeviation);

xAggregationScore = 1./sum(diff([0,find(xMeanDeviation < mean(xMeanDeviation))])>1);
yAggregationScore = 1./sum(diff([0;((find(yMeanDeviation < mean(yMeanDeviation))))])>1);


%%%%%% Plot Generation
subplot(1,3,1)
hold off
imagesc(Uint8ImageData) 
hold on
plot3(1:numel(xMeanDeviation),zeros(size(xMeanDeviation)),xMeanDeviation,'LineWidth',2)
plot3(zeros(size(yMeanDeviation))+numel(xMeanDeviation),1:numel(yMeanDeviation),yMeanDeviation,'LineWidth',2)
xlabel('X')
ylabel('Y')
axis([0,numel(xMeanDeviation),0,numel(yMeanDeviation),0,max(max(yMeanDeviation),max(xMeanDeviation))])
view(3)

subplot(1,3,2)
hold off
plot(xMeanDeviation)
hold on
plot(xMeanDeviation.*0+mean(xMeanDeviation))
xlabel('X')
ylabel('Relative Intensity')
title(sprintf('X-Aggregation Score = %g',xAggregationScore))
axis([0,numel(xMeanDeviation),0,max(xMeanDeviation)])

subplot(1,3,3)
hold off
plot(yMeanDeviation)
hold on
plot(yMeanDeviation.*0+mean(yMeanDeviation))
xlabel('Y')
ylabel('Relative Intensity')
title(sprintf('Y-Aggregation Score = %g',yAggregationScore))
axis([0,numel(yMeanDeviation),0,max(yMeanDeviation)])

end