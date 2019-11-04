function O = IspkIQs(IspkIbySpkCnt,TrlRange,QV)
O = quantile(IspkIbySpkCnt(TrlRange,:),QV);